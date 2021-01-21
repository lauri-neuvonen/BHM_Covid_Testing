#### Collection of tools to help with
# Building epidemic and optimization models
# Analysis of optimization and clustering results
# Most functions rely on core topical functions in
# ... policy_epidemic_model_code.py
# ... kMedoids_clustering.py


import numpy as np
import pandas as pd
from policy_epidemic_model_code import optimizable_corona_model

from run_definitions import *
from kMedoids_clustering import GetSwitch, GetCluster, CalcPearson, GetInitialMedoids, kMedoids

class Policy_template():

    def __init__(self, lockdown_policy_control_days, lockdown_policy_lower_limits,
                 lockdown_policy_upper_limits, testing_policy_control_days, testing_policy_lower_limits,
                 testing_policy_upper_limits):
        self.lockdown_policy_control_days = lockdown_policy_control_days
        self.lockdown_policy_lower_limits = lockdown_policy_lower_limits
        self.lockdown_policy_upper_limits = lockdown_policy_upper_limits
        self.testing_policy_control_days = testing_policy_control_days
        self.testing_policy_lower_limits = testing_policy_lower_limits
        self.testing_policy_upper_limits = testing_policy_upper_limits


class Policy():
    def __init__(self, lockdown_policy, testing_policy):
        self.lockdown_policy = lockdown_policy
        self.testing_policy = testing_policy

# Run generator

# NOTE: default values for all adjustable run parameters defined in function definition below:
def create_epidemic_model(ksi_base=ksi_base_default,
               A_rel=A_rel_default,
               r_AP=r_AP_default,
               r_U=r_U_default,
               r_P=r_P_default,
               r_N=r_N_default,
               d_vaccine=d_vaccine_default, # this is in _time steps_, not in days (days would be better)
               rel_rho=rel_rho_default,
               delta_param=delta_param_default,
               omegaR_param=omegaR_param_default,
               pii_D=pii_D_default,
               R_0=R_0_default,
               lambda_param = lambda_param_default,
               rel_lambda_param=rel_lambda_param_default,
               gamma_param=gamma_param_default,
               initial_infect=initial_infect_default,
               daily_testing_rate=daily_testing_rate_default,
               testing_sensitivity=testing_sensitivity_default,
               testing_specificity=testing_specificity_default,
               tau_TT_daily=tau_TT_daily_default,
               eta=eta_default,
               unknown_q_rate=unknown_q_rate_default,
               recovered_q_rate=recovered_q_rate_default,
               negative_q_rate=negative_q_rate_default,
               positive_q_rate=positive_q_rate_default,
               testing_cost=testing_cost_default,
               **kwargs):

    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, lambda_param, rel_lambda_param, initial_infect, testing_cost, eta, gamma_param)

    model_case = {
        'tau_paramA': (1 + daily_testing_rate) ** (1. / model.Delta_time) - 1,
        'test_sens': testing_sensitivity,
        'test_spec': testing_specificity,
        'tau_TT': (1 + tau_TT_daily) ** (1. / model.Delta_time) - 1,
        'ksi_U': (1 + unknown_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_P': (1 + positive_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_N': (1 + negative_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_R': (1 + recovered_q_rate) ** (1. / model.Delta_time) - 1,
        'r_U': (1 + r_U) ** (1. / model.Delta_time) - 1,
        'r_P': (1 + r_P) ** (1. / model.Delta_time) - 1,
        'r_AP': 0.0,
        'r_N': (1 + r_N) ** (1. / model.Delta_time) - 1,
        'r_R': (1 + r_N) ** (1. / model.Delta_time) - 1,  # same rate as for negatives
        'd_start_exp': 0.,
        'experiment': "baseline_vaccine_tag"
    }
    return model, model_case


def create_policy_control(lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    lockdown_policy_lower_limits=list(0.5 * np.ones(15)),
                    # can't use len(l_p_c_d) within function param def
                    lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
                    testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    testing_policy_lower_limits=list(np.zeros(15)),
                    testing_policy_upper_limits=list(0.02 * np.ones(15)),
                    ):

    policy_control = Policy_template(lockdown_policy_control_days, lockdown_policy_lower_limits,
                                     lockdown_policy_upper_limits, testing_policy_control_days,
                                     testing_policy_lower_limits,
                                     testing_policy_upper_limits)

    return policy_control


# NOTE: default values for all adjustable run parameters defined in function definition below:
def create_simu_run(ksi_base=ksi_base_default,
               A_rel=A_rel_default,
               r_AP=r_AP_default,
               r_U=r_U_default,
               r_P=r_P_default,
               r_N=r_N_default,
               d_vaccine=d_vaccine_default, # this is in _time steps_, not in days (days would be better)
               rel_rho=rel_rho_default,
               delta_param=delta_param_default,
               omegaR_param=omegaR_param_default,
               pii_D=pii_D_default,
               R_0=R_0_default,
                    lambda_param = lambda_param_default,
               rel_lambda_param=rel_lambda_param_default,
               gamma_param=gamma_param_default,
               initial_infect=initial_infect_default,
               daily_testing_rate=daily_testing_rate_default,
               testing_sensitivity=testing_sensitivity_default,
               testing_specificity=testing_specificity_default,
               tau_TT_daily=tau_TT_daily_default,
               eta=eta_default,
               unknown_q_rate=unknown_q_rate_default,
               recovered_q_rate=recovered_q_rate_default,
               negative_q_rate=negative_q_rate_default,
               positive_q_rate=positive_q_rate_default,
               testing_cost=testing_cost_default,


                lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                lockdown_policy_lower_limits=list(0.5 * np.ones(15)),
                # can't use len(l_p_c_d) within function param def
                lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
                testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                testing_policy_lower_limits=list(np.zeros(15)),
                testing_policy_upper_limits=list(0.02 * np.ones(15)),
                    **kwargs

                    ):
    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, lambda_param, rel_lambda_param, initial_infect, testing_cost, eta,
                                     gamma_param)

    model_case = {
        'tau_paramA': (1 + daily_testing_rate) ** (1. / model.Delta_time) - 1,
        'test_sens': testing_sensitivity,
        'test_spec': testing_specificity,
        'tau_TT': (1 + tau_TT_daily) ** (1. / model.Delta_time) - 1,
        'ksi_U': (1 + unknown_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_P': (1 + positive_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_N': (1 + negative_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_R': (1 + recovered_q_rate) ** (1. / model.Delta_time) - 1,
        'r_U': (1 + r_U) ** (1. / model.Delta_time) - 1,
        'r_P': (1 + r_P) ** (1. / model.Delta_time) - 1,
        'r_AP': 0.0,
        'r_N': (1 + r_N) ** (1. / model.Delta_time) - 1,
        'r_R': (1 + r_N) ** (1. / model.Delta_time) - 1,  # same rate as for negatives
        'd_start_exp': 0.,
        'experiment': "baseline_vaccine_tag"
    }

    policy_control = Policy_template(lockdown_policy_control_days, lockdown_policy_lower_limits,
                                     lockdown_policy_upper_limits, testing_policy_control_days,
                                     testing_policy_lower_limits,
                                     testing_policy_upper_limits)

    return model, model_case, policy_control


def create_sub_policy(policy_control_times, policy_control_values):
    policy = {}  # this will hold the policy in format suitable for input to the epidemic model
    #print("sub_policy: ", policy_control_times, " ", policy_control_values)
    if policy_control_times == "NA":
        return "NA"
    else:
        for (i, t) in enumerate(policy_control_times):
            policy[t] = policy_control_values[i]

        return policy

def create_policy(lockdown_policy, testing_policy):

    return Policy(lockdown_policy, testing_policy)

def cluster_run(run, run_policies_df, n_clusters):
    #run_control_times = list(map(int, run_policies_df.columns))
    run_policies = run_policies_df.to_numpy()

    run_obj_df = pd.read_csv('active_results/' + run + '_objectives.csv', delimiter=',')
    run_obj = run_obj_df.to_numpy()

    corr, dist = CalcPearson(run_policies)
    #print("dist: ", dist)
    cluster, medoids, cost = kMedoids(n_clusters, dist)
    return cluster, medoids, cost

def collect_results(runs, save_csv=False):
    # this function just combines solution and objective information from _results.csv and _objectives.csv files
    # into one dataframe and saves it as  full_results.csv if save_csv=True
    full_results_all = {}
    for run in runs:
        res_df = pd.read_csv('active_results/' + run + '_results.csv', delimiter=',')
        obj_df = pd.read_csv('active_results/' + run + '_objectives.csv', delimiter=',')
        full_results = res_df.join(obj_df)

        if save_csv:
            full_results.to_csv('active_results/' + run + '_full_results' + '.csv')

        full_results_all[run] = full_results

    return full_results_all

def extract_selected(runs, selected, save_csv=False, csv_identifier='selected'):

    # extracts selected rows from result dataframes corresponding to runs based on 'selected' list of dataframe indices
    # if save_csv is True a new csv is saved with identifier corresponding to csv_identifier. Full file name as below.
    selected_solutions = {}
    for run in runs:
        res_df = pd.read_csv('active_results/' + run + '_full_results.csv', delimiter=',', index_col=0)
        selected_res = res_df.loc[selected[run]]
        selected_solutions[run] = selected_res

        if save_csv:
            selected_res.to_csv('active_results/' + run + '_' + csv_identifier + '.csv')

    return selected_solutions


def simulate_solutions(run_list, result_set='full_results', no_control_runs=[]):
    policies = {}  # library to hold policy values, organized by: run, policy_number, policy values
    policy_obj_values = {}  # library to hold objective values, organized by: run, policy_number, objective values
    policy_sim_data = {}  # library to hold simulation output, organized by: run, policy_number, output_id, output_values
    no_control_sim_data = {}  # same as policy_sim_data except for there being no (different) policies
    epidemic_simulators = {}
    policy_controls = {}
    no_control_policy = Policy({10000: 0}, {10000: 0})

    for run in run_list:
        print("simulating policies for ", run)

        model, model_case, policy_control = create_simu_run(**run_list[run])

        epidemic_simulators[run] = (model, model_case)
        policy_controls[run] = policy_control  # ToDo: check if redundant
        # debug:
        # print("run: ", run)

        if run in no_control_runs:
            no_control_sim_data[run] = epidemic_simulators[run][0].solve_case(epidemic_simulators[run][1],
                                                                              no_control_policy)
        else:
            run_result_df = pd.read_csv('active_results/' + run + '_' + result_set + '.csv', delimiter=',', index_col=0)
            run_obj_df = run_result_df[['Deaths', 'Economic impact']]
            run_policies_df = run_result_df.drop(columns=['Deaths', 'Economic impact'])
            run_obj = run_obj_df.to_numpy()
            run_policies = run_policies_df.to_numpy()

            try:
                lockdown_only = (run_list[run]['testing_policy_control_days'] == "NA")
            except:
                lockdown_only = False  # if testing not defined in the run, default used -> testing used

            try:
                testing_only = (run_list[run]['lockdown_policy_control_days'] == "NA")
            except:
                testing_only = False

            if lockdown_only:
                ld_control_times = list(
                    map(int, [eval(tup)[1] for tup in run_policies_df.columns if (eval(tup)[0] == 'ld')]))

                policies[run] = {}
                policy_sim_data[run] = {}

                for pol_no, pol in enumerate(run_policies):
                    ld_pol = create_sub_policy(ld_control_times, pol)
                    test_pol = "NA"
                    policies[run][pol_no] = create_policy(ld_pol, test_pol)
                    policy_sim_data[run][pol_no] = epidemic_simulators[run][0].solve_case(epidemic_simulators[run][1],
                                                                                          policies[run])
            elif testing_only:

                test_control_times = list(
                    map(int, [eval(tup)[1] for tup in run_policies_df.columns if (eval(tup)[0] == 'test')]))
                policies[run] = {}
                policy_sim_data[run] = {}

                for pol_no, pol in enumerate(run_policies):
                    test_pol = create_sub_policy(test_control_times, pol)
                    ld_pol = "NA"
                    policies[run][pol_no] = create_policy(ld_pol, test_pol)
                    policy_sim_data[run][pol_no] = epidemic_simulators[run][0].solve_case(epidemic_simulators[run][1],
                                                                                          policies[run][pol_no])

            else:
                ld_control_times = list(
                    map(int, [eval(tup)[1] for tup in run_policies_df.columns if (eval(tup)[0] == 'ld')]))
                test_control_times = list(
                    map(int, [eval(tup)[1] for tup in run_policies_df.columns if (eval(tup)[0] == 'test')]))

                policies[run] = {}
                policy_sim_data[run] = {}

                for pol_no, pol in enumerate(run_policies):
                    ld_pol = create_sub_policy(ld_control_times, pol[0:len(ld_control_times)])
                    test_pol = create_sub_policy(test_control_times, pol[len(ld_control_times):len(pol)])

                    policies[run][pol_no] = create_policy(ld_pol, test_pol)
                    policy_sim_data[run][pol_no] = epidemic_simulators[run][0].solve_case(epidemic_simulators[run][1],
                                                                                          policies[run])

    print("Done.")
    return policies, policy_obj_values, policy_sim_data, epidemic_simulators