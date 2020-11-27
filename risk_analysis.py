import policy_epidemic_model_code
import numpy as np
import pandas as pd

from run_tools import create_epidemic_model, create_sub_policy, create_policy, Policy, Policy_template
from run_definitions import get_runs_definitions
from math import floor
import argparse
from run_definitions import (ksi_base_default, A_rel_default, r_AP_default, r_U_default, r_P_default, r_N_default, d_vaccine_default, rel_rho_default, delta_param_default,
    omegaR_param_default, pii_D_default, R_0_default, rel_lambda_param_default, gamma_param_default, initial_infect_default, daily_testing_rate_default,
    testing_sensitivity_default, testing_specificity_default, tau_TT_daily_default, eta_default, unknown_q_rate_default, recovered_q_rate_default,
    negative_q_rate_default, positive_q_rate_default, testing_cost_default)

# command line inputs
# first argument = sample size, integer
# second argument = runs to be analyzed, run names as strings, e.g. 'romer' 'romer_R0_4.0' ...
parser = argparse.ArgumentParser(description='Run optimization run using an NSGA-II algorithm for COVID-19 simulator')
parser.add_argument('sample_size', type=int, help='sample size for risk analysis')
parser.add_argument('runs', type=str, nargs='+', help='dictionary of run values to be changed from defaults')

args = parser.parse_args()
sample_size = args.sample_size # sample size from input
run_list = args.runs # runs to be analysed

run_definitions = get_runs_definitions()

runs = { run: run_definitions[run] for run in run_list } # filters the correct run definitions based on given run list


def sample_CVaR(sample, alpha, lowest_alpha=True):

    # this function calculates a CVaR -type 'worst case' expectation value for event mass corresponding to alpha
    # e.g. of sample group of N=100, this selects 10 worst values and calculates their average
    # lowest_alpha controls if the alpha share is selected from the low end (or the high end) of sample. If true the
    # lowest alpha_N values are used.

    sample_size = len(sample)
    print("sample size:", sample_size)
    print("alpha*sample size: ", sample_size*alpha, " --> ", floor(sample_size*alpha))
    alpha_N = floor(sample_size*alpha) # number of cases corresponding to share alpha

    sample.sort(reverse=lowest_alpha)
    worst_alpha = sample[-alpha_N:]

    remainder = (sample_size*alpha - alpha_N)*sample[-(alpha_N+1)]
    print("worst alpha: ", worst_alpha, " remainder: ", remainder, 'from', sample[-(alpha_N+1)])

    return (sum(worst_alpha)+remainder)/(sample_size*alpha)

for run in runs:

    # get policies:

    run_policies_df = pd.read_csv('active_results/' + run + '_results.csv', delimiter=',')
    run_control_times = list(map(int, run_policies_df.columns))
    run_policies = run_policies_df.to_numpy()


    # get parameter value samples:
    print("run: ", runs[run])
    sample_list = []

    ### INPUT RISK ANALYSIS DEFINITIONS BELOW! ###
    ### MAKE SURE CORRECT SOURCE FOR EXP is used ###

    analysis_params = {}

    try:
        analysis_params['R_0']= {
            '_dist': np.random.lognormal,
            'dist_params': (np.log(runs[run]['R_0']**2/np.sqrt((runs[run]['R_0']**2 + 0.2))), np.log(1 + 0.2/runs[run]['R_0']**2))  # [mean, sigma] for lognormal
        }
    except KeyError:
        analysis_params['R_0'] = {
            '_dist': np.random.lognormal,
            'dist_params': (np.log(R_0_default**2/np.sqrt((R_0_default**2 + 0.2))), np.log(1 + 0.2/R_0_default**2))   # [mean, sigma] for lognormal
        }

    analysis_params['pii_D'] = {
        '_dist': np.random.lognormal,
        'dist_params': (
        np.log(pii_D_default ** 2 / np.sqrt((pii_D_default ** 2 + 0.0002))), np.log(1 + 0.0002 / pii_D_default ** 2))
        # [mean, sigma] for lognormal
    }

    for i in range(0, sample_size):
        sample_instance = {}
        for param in analysis_params:
            sample_instance[param] = analysis_params[param]['_dist'](*analysis_params[param]['dist_params'])
        sample_list.append(sample_instance)
   #debug
    print("samples: ", sample_list)


    run_policy_samples = {}

    for policy_id, policy in enumerate(run_policies):

        # create policy for run:

        try:
            lockdown_only = (runs[run]['testing_policy_control_days'] == "NA")
        except:
            lockdown_only = False  # if testing not defined in the run, default used -> testing used

        try:
            testing_only = (runs[run]['lockdown_policy_control_days'] == "NA")
        except:
            testing_only = False

        if lockdown_only:
            ld_policy = create_sub_policy(run_control_times, policy)
            test_policy = "NA"

        elif testing_only:
            test_policy = create_sub_policy(run_control_times, policy)
            ld_policy = "NA"

        else:
            ld_control_times = run_control_times[:len(run_control_times)//2]
            ld_policy = create_sub_policy(ld_control_times, policy[:len(ld_control_times)])

            test_control_times = run_control_times[len(run_control_times)//2:]
            test_policy = create_sub_policy(test_control_times, policy[len(ld_control_times):])

        run_policy = create_policy(ld_policy, test_policy)

        policy_result_dist = {}
        policy_sample_deaths = []

        for sample_id, sample in enumerate(sample_list):
            sample_run_params = runs[run].copy() # copies the original run (e.g. 'romer') for updating with sample values
            sample_run_params.update(sample)
            #debug:
            #print("sample params: ", sample_run_params)

            # calculate results

            epidemic_simulator = create_epidemic_model(**sample_run_params)
            Reported_D, Notinfected_D, Unreported_D, Infected_D, \
            False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total, total_testing_cost, tests, Unk_NA_nQ_D, Unk_NA_Q_D, K_NA_nQ_D, Unk_IA_nQ_D, Unk_IA_Q_D, K_IA_Q_D, alpha_D, ksi_TT_I_D, ksi_TT_N_D, ksi_TT_R_D, Symptomatic_T \
                = epidemic_simulator[0].solve_case(epidemic_simulator[1], run_policy)

            #debug
            #print("simu deaths: ", Dead_D[-1], " | simu output: ", Y_total)

            policy_result_dist[sample_id] = [Dead_D[-1], Y_total]
            policy_sample_deaths.append(Dead_D[-1])

        print("deaths in sample: ", policy_sample_deaths)
        print("policy CVaR 10%: ", 340000000*sample_CVaR(policy_sample_deaths, 0.1, False))

        #run_policy_samples[sample_id] = policy_result_dist

        # saving results to csv files
        df = pd.DataFrame(policy_result_dist, index=['Deaths', 'Output'])

        df.to_csv('active_results/risk_analysis/'+run+'__'+str(policy_id)+'.csv')



