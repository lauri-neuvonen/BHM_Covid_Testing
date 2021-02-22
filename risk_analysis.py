import numpy as np
from scipy.stats import truncnorm
import pandas as pd

from run_tools import create_epidemic_model, create_sub_policy, create_policy
from run_definitions import get_runs_definitions, p_ICU_def, C_hos_def, T_rec_def
from math import floor
import argparse
from progress.bar import Bar

from run_definitions import initial_infect_default

# command line inputs
# first argument = sample size, integer
# second argument = runs to be analyzed, run names as strings, e.g. 'romer' 'romer_R0_4.0' ...
parser = argparse.ArgumentParser(description='Risk analysis using random generated parameter samples.')
parser.add_argument('sample_size', type=int, help='sample size for risk analysis')
parser.add_argument('result_set', type=str, help='selects which result file to use: full_results or medoid')
parser.add_argument('runs', type=str, nargs='+', help='dictionary of run values to be changed from defaults')
parser.add_argument('--file_suffix', type=str, help='suffix to add to the end on filename')
parser.add_argument('--policy_file', type=str, help='Optional. If set, policies are read from this file instead of one determined by run and result set id.')
parser.add_argument('--params', type=str, nargs='+', help='parameters to include in sensitivity analysis, see definitions for all available.')
#parser.add_argument('--policy_index', type=int, help='Optional. If set, the corresponding policy from policy_file is used as policy for sample simulations')


args = parser.parse_args()
sample_size = args.sample_size # sample size from input
run_list = args.runs # runs to be analysed
set_id = args.result_set
file_suffix = args.file_suffix
policy_file = args.policy_file
params_sel = args.params # list of parameters to include in sensitivity analysis

run_definitions = get_runs_definitions()

runs = { run: run_definitions[run] for run in run_list } # filters the correct run definitions based on given run list

def trunc_norm_builder(my_lower, my_upper, mu, sigma):
    lower  = (my_lower - mu) / sigma
    upper = (my_upper - mu) / sigma
    return truncnorm(lower, upper, loc=mu, scale=sigma).rvs(1)[0]

def sample_CVaR(sample, alpha, lowest_alpha=True):

    # this function calculates a CVaR -type 'worst case' expectation value for event mass corresponding to alpha
    # e.g. of sample group of N=100, this selects 10 worst values and calculates their average
    # lowest_alpha controls if the alpha share is selected from the low end (or the high end) of sample. If true the
    # lowest alpha_N values are used.

    sample_size = len(sample)
    #print("sample size:", sample_size)
    #print("alpha*sample size: ", sample_size*alpha, " --> ", floor(sample_size*alpha))
    alpha_N = floor(sample_size*alpha) # number of cases corresponding to share alpha

    sample.sort(reverse=lowest_alpha)
    if alpha_N !=0:    # makes sure -0 doesn't cause trouble in calculating worst_alpha
        worst_alpha = sample[-alpha_N:]
    else:
        worst_alpha = []

    remainder = (sample_size*alpha - alpha_N)*sample[-(alpha_N+1)]
    #print("worst alpha: ", worst_alpha, " remainder: ", remainder, 'from', sample[-(alpha_N+1)])

    return (sum(worst_alpha)+remainder)/(sample_size*alpha)


# get parameter value samples:
#print("run: ", runs[run])
sample_list = []

### INPUT RISK ANALYSIS DEFINITIONS BELOW! ###

analysis_params = {}
rng = np.random.default_rng(12345) # creates a random number generator with constant seed sequence

analysis_params['R_0']= {
    '_dist': rng.gamma,
    'dist_params': (100,0.025)  # (low, high) for uniform
}

analysis_params['pii_D'] = {
    '_dist': rng.beta,
    'dist_params': (1.45, 95) # (low, high) for uniform
}

analysis_params['delta_param']= {
    '_dist': rng.gamma,
    'dist_params': (2.81, 2.385)
}

analysis_params['gamma_param'] = {
    '_dist': rng.uniform,
    'dist_params': (150, 240)
}

analysis_params['initial_infect'] = {
    '_dist': rng.uniform,
    'dist_params': (initial_infect_default / 2, initial_infect_default * 1.5)
}


if params_sel == None:  # if no params are listed, default is to use all.
    params_sel = analysis_params.keys()


param_values = {}


for p in params_sel:  # creates an empty array for each param for storing sample values
    param_values[p] = []

for i in range(0, sample_size):
    sample_instance = {}
    for param in params_sel:
        sample_instance[param] = analysis_params[param]['_dist'](*analysis_params[param]['dist_params']) # saves param value to sample point
        param_values[param].append(sample_instance[param]) # saves the param value to param specific list

    sample_list.append(sample_instance)

for run in runs:

    if policy_file != None:
        file_path = policy_file
    else:
        file_path = 'active_results/' + run + "_" + set_id + '.csv'

    print("Downloads policies from: ", file_path)
        # get policies:

    # control times saved as tuples, e.g. ('ld', 100) for 'lockdown at time 100'.
    run_results_df = pd.read_csv(file_path, delimiter=',')
    run_policies_df = run_results_df.drop(columns=['Deaths', 'Economic impact', 'cluster'], errors='ignore')

    # next lines extract control times from column names by evaluating the string into a tuple and picking out
    # the times from the tuple. The times are then mapped to integers and listed.
    ld_control_times = list(map(int, [eval(tup)[1] for tup in run_policies_df.columns if eval(tup)[0] == 'ld']))
    test_control_times = list(map(int, [eval(tup)[1] for tup in run_policies_df.columns if eval(tup)[0] == 'test']))

    run_policy_samples = {}
    policy_CVaRs = []
    policy_ICUOL_Ps = []


    bar = Bar('Simulating policies', max=len(run_policies_df.index))
    for policy_id, policy in run_policies_df.iterrows():
        #print("\npolicy id: ", policy_id, ": ", policy) # debug

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
            ld_policy = create_sub_policy(ld_control_times, policy)
            test_policy = "NA"

        elif testing_only:
            test_policy = create_sub_policy(test_control_times, policy)
            ld_policy = "NA"

        else:
            # ld_control_times = run_control_times[:len(run_control_times)//2]
            ld_policy = create_sub_policy(ld_control_times, policy[:len(ld_control_times)])

            # test_control_times = run_control_times[len(run_control_times)//2:]
            test_policy = create_sub_policy(test_control_times, policy[len(ld_control_times):])


        run_policy = create_policy(ld_policy, test_policy)

        policy_result_dist = {}
        policy_sample_ICUover = []
        policy_sample_ICU_bool = []

        for sample_id, sample in enumerate(sample_list):

            sample_run_params = runs[run].copy() # copies the original run (e.g. 'romer') for updating with sample values
            sample_run_params.update(sample)

            # calculate results

            epidemic_simulator = create_epidemic_model(**sample_run_params)
            Reported_D, Notinfected_D, Unreported_D, Infected_D, \
            False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total, total_testing_cost, tests, Unk_NA_nQ_D, Unk_NA_Q_D, K_NA_nQ_D, Unk_IA_nQ_D, Unk_IA_Q_D, K_IA_Q_D, alpha_D, ksi_TT_I_D, ksi_TT_N_D, ksi_TT_R_D, Symptomatic_T \
                = epidemic_simulator[0].solve_case(epidemic_simulator[1], run_policy)

            # calculating aggregated ICU capacity overload:

            try:
                p_ICU = runs[run]['p_ICU']
            except:
                p_ICU = p_ICU_def

            try:
                C_hos= runs[run]['C_hos']
            except:
                C_hos = C_hos_def

            ICU_use_T = p_ICU * Symptomatic_T
            ICU_margin_T = ICU_use_T - (C_hos / epidemic_simulator[0].pop)
            ICU_overuse_T = np.max([ICU_margin_T], initial=0.0, axis=0) # this creates a vector which has non-zero values for overuse, 0 otherw.
            ICU_overuse_agg = np.sum(ICU_overuse_T)


            try:
                T_rec = runs[run]['T_rec']
            except:
                T_rec = T_rec_def

            T_rec_t = int(round(14 * 365 * T_rec))  # change from years to time steps
            # Costs:
            cost_e = -Y_total / epidemic_simulator[0].T  # contains loss of output & scaled direct costs
            cost_terminal = ((T_rec_t) / 2) * (-Y_D[-1]) / epidemic_simulator[0].T
            deaths_terminal = ((T_rec_t) / 2) * ((Dead_D[-1] - Dead_D[
                -2]) * epidemic_simulator[0].pop / 1000) / epidemic_simulator[0].T  # Deaths are cumulative, so difference needed for current rate
            # hcap_terminal = ((T_rec_t) / 2) * np.max([0.0, self.p_ICU * Symptomatic_T[-1] - self.C_hos / self.model.pop])
            # print("cost_e: ", cost_e)
            # print("cost_terminal: ", cost_terminal)

            # objectives scaled in same way as in optimizer objective calculation!
            scaling = epidemic_simulator[0].T_years / (T_rec + epidemic_simulator[0].T_years)
            deaths_norm_adj = Dead_D[-1] * epidemic_simulator[0].pop / 1000 + deaths_terminal
            output_norm_adj = scaling * (cost_e + cost_terminal)

            policy_result_dist[sample_id] = [deaths_norm_adj, output_norm_adj, ICU_overuse_agg]
            policy_sample_ICUover.append(ICU_overuse_agg)

            ICU_overuse_T_bool = list(map(bool, ICU_overuse_T))
            policy_sample_ICU_bool.append(ICU_overuse_T_bool)

        #print("overload in sample: ", policy_sample_ICUover)
        policy_CVaR = 340000000*sample_CVaR(policy_sample_ICUover, 0.1, False)
        policy_CVaRs.append(policy_CVaR)
        #print("policy CVaR 10%: ", 340000000*sample_CVaR(policy_sample_ICUover, 0.1, False))

        #print("policy sample ICU bool: ", policy_sample_ICU_bool)
        max_sample_overloads_T = max(np.sum(policy_sample_ICU_bool, axis=0))
        max_P_ICU_overload = max_sample_overloads_T / sample_size
        policy_ICUOL_Ps.append(max_P_ICU_overload)

        # saving results to csv files
        df = pd.DataFrame.from_dict(policy_result_dist, orient='index', columns=['Deaths', 'Output', 'ICU overload'])
        #print("df: ", df)
        #print("sample: ", sample_pii_Ds)
        for p in params_sel: # inserts the param value lists to df
            df.insert(0, p, param_values[p])


        df.to_csv('active_results/risk_analysis/'+run+'__'+str(policy_id) + set_id + '_' + file_suffix + '.csv')
        bar.next()

    # get full results:

    full_results = run_results_df.copy()

    loc = len(full_results.columns)
    full_results.insert(loc, 'ICUOL (CVaR 10%)', policy_CVaRs)
    full_results.insert(loc+1, 'max ICUOL P', policy_ICUOL_Ps)
    full_results.to_csv('active_results/risk_analysis/'+run+'_' + set_id + '_risk' + file_suffix +'.csv')

