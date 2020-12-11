from policy_epidemic_model_code import *
from jupyterWidgets import *

import numpy as np
import importlib
import matplotlib.pyplot as plt

import pandas as pd

#
#   , alpha_D, ksi_TT_I_D, ksi_TT_N_D, ksi_TT_R_D, Symptomatic_T
output_names = {
    0: "Reported_D",
    1: "Notinfected_D",
    2: "Unreported_D",
    3: "Infected_D",
    4: "False_pos (end of day)",
    5: "False_neg (end of day)",
    6: "Recovered_D",
    7: "Dead_D",
    8: "Infected_T (time step)",
    9: "Infected_not_Q",
    10: "Infected_in_Q",
    11: "Output, (end of day)",
    12: "M_t",
    13: "Y_total",
    14: "total_cost",
    15: "tests",
    16: "Unknown, not infected, asymptomatic, not quarantined",
    17: "Unknown, not infected, asymptomatic, quarantined",
    18: "Known, not infected, asymptomatic, not quarantined",
    19: "Unknown, infected, asymptomatic, not quarantined",
    20: "Unknown, infected, asymptomatic, quarantined",
    21: "Known, infected, asymptomatic, quarantined",
    22: "alpha_T",
    23: "ksi_TT_I_T",
    24: "ksi_TT_N_D",
    25: "ksi_TT_R_D",
    26: "Symptomatic_T"
}

run_labels = {
    'base_case_no_control':         'no control',
    'base_case_no_control_R0_4.0':      'no control with R0=4.0',
    'romer':        'mass testing',
    'romer_no_limit':       'mass testing no testing limit',
    'romer_terminal_0.25':      'mass testing, 0.25y terminal value',
    'romer_terminal_1.0':       'mass testing, 1y terminal value',
    'romer_R0_4.0':     'mass testing, R0=4.0',
    'romer_R0_4.0_no_limit':        'mass testing, R0=4.0, no testing limit',
    'romer_R0_1.25':        'mass testing, R0=1.25',
    'romer_R0_4.0_sens_spec_075':       'mass testing, R0=4.0, sens&spec 0.75',
    'romer_R0_4.0_sens_spec_085':       'mass testing, R0=4.0, sens&spec 0.85',
    'romer_R0_4.0_sens_spec_090':       'mass testing, R0=4.0, sens&spec 0.90',
    'romer_R0_4.0_sens_spec_095':       'mass testing, R0=4.0, sens&spec 0.95',
    'romer_6d_incubation':      'mass testing, 6d incub.',
    'romer_6d_incubation_sens_spec_090':        'mass testing, 6d incub., sens&spec 0.90',
    'romer_8d_incubation':      'mass testing, 8d incub.',
    'romer_8d_incubation_sens_spec_075':        'mass testing, 8d incub., sens&spec 0.75',
    'romer_8d_incubation_sens_spec_090':        'mass testing, 8d incub., sens&spec 0.90',
    'romer_3d_delay':       'mass testing, 3d delay',
    'romer_7d_delay':       'mass testing, 7d delay',
    'romer_14d_delay':      'mass testing, 14d delay',
    'romer_28d_delay':      'mass testing, 28d delay',
    'romer_sens_085':       'mass testing, sens 0.85',
    'romer_spec_085':       'mass testing, spec 0.85',
    'romer_sens_spec_085':      'mass testing, sens&spec 0.85',
    'romer_sens_spec_090':      'mass testing, sens&spec 0.90',
    'romer_sens_spec_095':      'mass testing, sens&spec 0.95',
    'romer_sens_spec_098':      'mass testing, sens&spec 0.98',
    'romer_sens_spec_099':      'mass testing, sens&spec 0.99',
    'base_case_lockdown_opt':   'lockdown',
    'base_case_lockdown_opt_terminal_0.25':     'lockdown, 0.25y terminal value',
    'base_case_lockdown_opt_terminal_1.0':     'lockdown, 1y terminal value',
    'base_case_lockdown_opt_R0_4.0':    'lockdown, R0=4.0',
    'base_case_lockdown_opt_with_limited_general_testing': 'lockdown, w. perf. testing',
    'base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing':  'lockdown, w. testing sens&spec=0.75',
    'base_case_lockdown_opt_with_limited_sens075_general_testing':  'lockdown, w. testing sens=0.75',
    'base_case_lockdown_opt_with_limited_spec075_general_testing':  'lockdown, w. testing spec=0.75',
    'base_case_lockdown_opt_with_limited_sens090_general_testing':  'lockdown, w. testing sens=0.90',
    'base_case_lockdown_opt_with_limited_spec090_general_testing':  'lockdown, w. testing spec=0.90',
    'base_case_lockdown_opt_with_limited_imperfect(0.85)_general_testing':  'lockdown, w. testing sens&spec=0.85',
    'base_case_lockdown_opt_with_limited_imperfect(0.90)_general_testing':  'lockdown, w. testing sens&spec=0.90',
    'base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing':  'lockdown, w. testing sens&spec=0.95',
    'base_case_lockdown_opt_3d_delay':  'lockdown, 3d delay',
    'base_case_lockdown_opt_7d_delay':  'lockdown, 7d delay',
    'base_case_lockdown_opt_14d_delay':  'lockdown, 14d delay',
    'base_case_lockdown_opt_28d_delay':  'lockdown, 28d delay',
    'base_case_6d_incubation': 'lockdown 6d incub.',
    'base_case_8d_incubation': 'lockdown 8d incub.',
    'test_and_trace_no_control_eta50': 'no control w CT, eta=50',
    'test_and_trace_no_control_eta100':'no control w CT, eta=100',
    'test_and_trace_no_control_R0_4.0_eta10': 'no control, R0=4.0, w CT, eta=10',
    'test_and_trace_no_control_R0_4.0_eta50': 'no control, R0=4.0, w CT, eta=50',
    'test_and_trace_no_control_R0_4.0_eta100': 'no control, R0=4.0, w CT, eta=100',
    'test_and_trace_lockdown_opt_eta10':    'lockdown w CT, eta=10',
    'test_and_trace_lockdown_opt_eta50':    'lockdown w CT, eta=50',
    'test_and_trace_lockdown_opt_eta75':    'lockdown w CT, eta=75',
    'test_and_trace_lockdown_opt_eta95':    'lockdown w CT, eta=95',
    'test_and_trace_lockdown_opt_eta100':    'lockdown w CT, eta=100',
    'test_and_trace_lockdown_opt_eta50_R04':    'lockdown w CT, R0=4.0, eta=50',
    'test_and_trace_lockdown_opt_eta75_R04':    'lockdown w CT, R0=4.0, eta=75',
    'test_and_trace_lockdown_opt_eta100_R04':    'lockdown w CT, R0=4.0, eta=100',
    'test_and_trace_lockdown_opt_eta50_R04_delta10':    'lockdown w CT, R0=4.0, eta=50, delta=10d',
    'test_and_trace_lockdown_opt_eta75_R04_delta10':    'lockdown w CT, R0=4.0, eta=75, delta=10d',
    'test_and_trace_lockdown_opt_eta100_R04_delta10':    'lockdown w CT, R0=4.0, eta=100, delta=10d',
    'combo_base_case': 'combo',
    'combo_base_case_R0_4.0': 'combo, R0=4.0',
    'combo_base_case_test_and_trace': 'combo CT',
    'combo_sens_spec_0.95': 'combo, sens&spec=0.95',
    'combo_sens_spec_0.85': 'combo, sens&spec=0.85',
    'combo_R0_4.0_sens_spec_0.95': 'combo, R0=4.0, sens&spec=0.95',
    'combo_R0_4.0_sens_spec_0.85': 'combo, R0=4.0, sens&spec=0.85',
    'base_case_no_control_1000': 'no control, 1000 initial infected',
    'base_case_no_control_2000': 'no control, 2000 initial infected',
    'base_case_no_control_10000': 'no control, 10000 initial infected'
}


def epidemic_progression_plot(outputs, epidemic_sims, runs_data, columns=2, policies="NA"):

    rows = int(np.ceil(len(outputs)/columns))
    #empty = rows*columns - len(outputs)
    #fills = [-1]*empty
    #outputs_filled = outputs.extend(fills)
    outputs_array = np.reshape(outputs, (rows, columns))
    fig, axes = plt.subplots(nrows=rows, ncols=columns, figsize=(18, 9))


    for row in range(0, rows):
        for col in range(0, columns):
            if policies != "NA":
                secax = axes[row, col].twinx()
            for run in runs_data:
                time_steps = range(0, epidemic_sims[run][0].T_years*365)
                axes[row,col].plot(time_steps, runs_data[run][outputs_array[row,col]], label=run_labels[run])
                axes[row,col].set_xlabel('time (days)')
                #axes[row,col].set_ylabel(output_names[outputs_array[row,col]])

                if policies != "NA":
                    try:
                        secax.step(list(policies[run].testing_policy.keys()), list(policies[run].testing_policy.values()), ":", where='post',
                                   label=run_labels[run] + ": testing policy")
                        secax.set_ylabel('policy variable value')
                    except:
                        pass
                    try:
                        secax.step(list(policies[run].lockdown_policy.keys()), list(policies[run].lockdown_policy.values()), ":", where='post',
                               label=run_labels[run] + ": lockd. policy")
                        secax.set_ylabel('policy variable value')
                    except:
                        pass

            axes[row, col].set_title(output_names[outputs_array[row,col]])
            # axes[row, col].legend()
            if policies != "NA":
                secax.legend()

    return fig

def pareto_plot(runs, path="active_results/"):

    #fig, axes = plt.subplots(ncols=2, figsize=(16, 8))
    fig, axes = plt.subplots(ncols=1, figsize=(8, 8))

    for run in runs:
        obj = pd.read_csv(path+run+"_objectives.csv", delimiter=',').to_numpy()
        axes.scatter(obj[:,0], -obj[:,1], label=run_labels[run])
        axes.set_title('Cumulative deaths vs total economic output')
        axes.set_xlabel('deaths, 1000 persons')
        axes.set_ylabel('total output (Y_total + terminal value)')
        axes.legend()

        #axes[1].scatter(obj[:, 0], obj[:, 2], label=run_labels[run])
        #axes[1].set_title('Cumulative deaths vs cost-output efficiency measure')
        #axes[1].set_xlabel('deaths, 1000 persons')
        #axes[1].set_ylabel('hospital capacity overload')
        #axes[1].legend()
    return fig
# Tools for building optimization runs based on params.

def create_policy(policy_control_times, policy_control_values):
    policy = {}  # this will hold the policy in format suitable for input to the epidemic model
    # print("times: ", policy_control_times)
    # print("values: ", policy_control_values)

    for (i, t) in enumerate(policy_control_times):
        policy[t] = policy_control_values[i]

    return policy


# Run generator

# NOTE: default values for all adjustable run parameters defined in function definition below:
def create_simu_run(ksi_base=0,
                    A_rel=0.5,
                    r_AP=0,
                    d_vaccine=800,
                    rel_rho=1.0,
                    delta_param=5,
                    omegaR_param=14,
                    gamma_param=180,
                    pii_D=0.01,
                    R_0=2.5,
                    rel_lambda_param=0.5,
                    initial_infect=300,
                    testing_rate=0.0,
                    testing_sensitivity=1.0,
                    testing_specificity=1.0,
                    tau_TT=0.0,
                    eta=0.0,
                    unknown_q_rate=0.0,
                    recovered_q_rate=0.0,
                    negative_q_rate=0.0,
                    positive_q_rate=0.999,
                    testing_cost=100,
                    pop_size=28,
                    n_offsprings=14,
                    # sampling=get_sampling("real_random"),
                    # crossover=get_crossover("real_sbx", prob=0.9, eta=15),
                    # mutation=get_mutation("real_pm", eta=15),
                    # eliminate_duplicates=True,
                    # filename = "foo",
                    # termination = get_termination("n_gen", 100),
                    lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    lockdown_policy_lower_limits=list(0.5 * np.ones(15)),
                    # can't use len(l_p_c_d) within function param def
                    lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
                    testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    testing_policy_lower_limits=list(np.zeros(15)),
                    testing_policy_upper_limits=list(0.2 * np.ones(15))
                    ):
    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect, testing_cost, eta,
                                     gamma_param)

    model_case = {
        'tau_paramA': testing_rate,
        'test_sens': testing_sensitivity,
        'test_spec': testing_specificity,
        'tau_TT': tau_TT,
        'ksi_U': (1 + unknown_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_P': (1 + positive_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_N': (1 + negative_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_R': (1 + recovered_q_rate) ** (1. / model.Delta_time) - 1,
        'r_U': (1 + 0.98) ** (1. / model.Delta_time) - 1,  # should be redundant!
        'r_P': (1 + 0.98) ** (1. / model.Delta_time) - 1,
        'r_AP': 0.0,
        'r_N': (1 + 0.98) ** (1. / model.Delta_time) - 1,
        'r_R': (1 + 0.999) ** (1. / model.Delta_time) - 1,
        'd_start_exp': 0.,
        'experiment': "baseline_vaccine_tag"
    }

    policy_control = Policy_template(lockdown_policy_control_days, lockdown_policy_lower_limits,
                                     lockdown_policy_upper_limits, testing_policy_control_days,
                                     testing_policy_lower_limits,
                                     testing_policy_upper_limits)

    return model, model_case, policy_control
