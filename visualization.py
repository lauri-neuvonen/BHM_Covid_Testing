from policy_epidemic_model_code import *
from jupyterWidgets import *

import numpy as np
import importlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import pandas as pd

#
#   , alpha_D, ksi_TT_I_D, ksi_TT_N_D, ksi_TT_R_D, Symptomatic_T
output_names = {
    0: "Reported_D",
    1: "Notinfected_D",
    2: "Unreported_D",
    3: "Infected_D",
    4: "False_pos",
    5: "False_neg",
    6: "Recovered_D",
    7: "Dead (cumulative)",
    8: "Infected_T (time step)",
    9: "Infected_not_Q",
    10: "Infected_in_Q",
    11: "Output",
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
    22: "Probability of infection from a random encounter",
    23: "ksi_TT_I_T",
    24: "ksi_TT_N_D",
    25: "ksi_TT_R_D",
    26: "Symptomatic"
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
    'combo_sens_spec_0.95': 'combo, sens&spec=0.95',
    'combo_sens_spec_0.85': 'combo, sens&spec=0.85',
    'combo_R0_4.0_sens_spec_0.95': 'combo, R0=4.0, sens&spec=0.95',
    'combo_R0_4.0_sens_spec_0.85': 'combo, R0=4.0, sens&spec=0.85',
    'base_case_no_control_1000': 'no control, 1000 initial infected',
    'base_case_no_control_2000': 'no control, 2000 initial infected',
    'base_case_no_control_10000': 'no control, 10000 initial infected',
    'romer_tc_1000000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 1M tests',
    'romer_tc_2500000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 2.5M tests',
    'romer_tc_5000000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 5M tests',
    'romer_tc_25000000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 25M tests',
    'romer_tc_50000000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 50M tests',
    'romer_tc_100000000_sens_spec_085': 'mass_testing, sens&spec=0.85, max 100M tests',
    'romer_tc_1000000': 'mass_testing, max 1M tests',
    'romer_tc_2500000': 'mass_testing, max 2.5M tests',
    'romer_tc_5000000': 'mass_testing, max 5M tests',
    'romer_tc_25000000': 'mass_testing, max 25M tests',
    'romer_tc_50000000': 'mass_testing, max 50M tests',
    'romer_tc_100000000': 'mass_testing, max 100M tests',
    'combo_base_case_test_and_trace': 'combo test&trace',
    'combo_base_case_tc_1000000': 'combo max 1M tests',
    'combo_base_case_tc_2500000': 'combo max 2.5M tests',
    'combo_base_case_tc_5000000': 'combo max 5M tests',
    'combo_base_case_tc_25000000': 'combo max 25M tests',
    'combo_base_case_tc_50000000': 'combo max 50M tests',
    'combo_base_case_tc_100000000': 'combo max 100M tests',
    'combo_base_case_test_and_trace_tc1000000': 'combo max 1M tests + test&trace',
    'combo_base_case_test_and_trace_tc2500000': 'combo max 2.5M tests + test&trace',
    'combo_base_case_test_and_trace_tc5000000': 'combo max 5M tests + test&trace',
    'combo_base_case_test_and_trace_tc25000000': 'combo max 25M tests + test&trace',
    'combo_base_case_test_and_trace_tc50000000': 'combo max 50M tests + test&trace',
    'combo_base_case_test_and_trace_tc100000000': 'combo max 100M tests + test&trace',
    'combo_base_case_test_and_trace_ss085': 'combo max 3M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc1000000_ss085': 'combo max 1M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc2500000_ss085': 'combo max 2.5M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc5000000_ss085': 'combo max 5M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc25000000_ss085': 'combo max 25M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc50000000_ss085': 'combo max 50M tests + test&trace, sens&spec=0.85',
    'combo_base_case_test_and_trace_tc100000000_ss085': 'combo max 100M tests + test&trace, sens&spec=0.85',
    'combo_base_case_tc_50000000_sens_spec_085': 'combo max 50M tests, sens&spec=0.85',
    'combo_base_case_tc_25000000_sens_spec_085': 'combo max 25M tests, sens&spec=0.85',
    'combo_base_case_tc_5000000_sens_spec_085': 'combo max 5M tests, sens&spec=0.85',
    'combo_base_case_tc_2500000_sens_spec_085': 'combo max 2.50M tests, sens&spec=0.85',
    'combo_base_case_tc_1000000_sens_spec_085': 'combo max 1M tests, sens&spec=0.85',
    'combo_base_case_tc_100000000_sens_spec_085': 'combo max 100M tests, sens&spec=0.85',
    'combo_base_case_sens_spec_085': 'combo max 3M tests, sens&spec=0.85'


}

short_run_labels = {
    'base_case_no_control':         'no control',
    'base_case_no_control_R0_4.0':      'no control with R0=4.0',
    'romer':        'MT, max 3M tests',
    'romer_no_limit':       'MT no testing limit',
    'romer_terminal_0.25':      'MT, max 3M tests, 0.25y terminal value',
    'romer_terminal_1.0':       'MT, max 3M tests, 1y terminal value',
    'romer_R0_4.0':     'MT, max 3M tests, R0=4.0',
    'romer_R0_4.0_no_limit':        'MT, max 3M tests, R0=4.0, no testing limit',
    'romer_R0_1.25':        'MT, max 3M tests, R0=1.25',
    'romer_R0_4.0_sens_spec_075':       'MT, R0=4.0, s&s 0.75',
    'romer_R0_4.0_sens_spec_085':       'MT, R0=4.0, s&s 0.85',
    'romer_R0_4.0_sens_spec_090':       'MT, R0=4.0, s&s 0.90',
    'romer_R0_4.0_sens_spec_095':       'MT, R0=4.0, s&s 0.95',
    'romer_6d_incubation':      'MT, 6d incub.',
    'romer_6d_incubation_sens_spec_090':        'MT, 6d incub., s&s 0.90',
    'romer_8d_incubation':      'MT, 8d incub.',
    'romer_8d_incubation_sens_spec_075':        'MT, 8d incub., s&s 0.75',
    'romer_8d_incubation_sens_spec_090':        'MT, 8d incub., s&s 0.90',
    'romer_3d_delay':       'MT, 3d delay',
    'romer_7d_delay':       'MT, 7d delay',
    'romer_14d_delay':      'MT, 14d delay',
    'romer_28d_delay':      'MT, 28d delay',
    'romer_sens_085':       'MT, sens 0.85',
    'romer_spec_085':       'MT, spec 0.85',
    'romer_sens_spec_085':      'MT, s&s 0.85',
    'romer_sens_spec_090':      'MT, s&s 0.90',
    'romer_sens_spec_095':      'MT, s&s 0.95',
    'romer_sens_spec_098':      'MT, s&s 0.98',
    'romer_sens_spec_099':      'MT, s&s 0.99',
    'base_case_lockdown_opt':   'LD, no testing',
    'base_case_lockdown_opt_terminal_0.25':     'LD, 0.25y terminal value',
    'base_case_lockdown_opt_terminal_1.0':     'LD, 1y terminal value',
    'base_case_lockdown_opt_R0_4.0':    'LD, R0=4.0',
    'base_case_lockdown_opt_with_limited_general_testing': 'LD, w. perf. testing',
    'base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing':  'LD, w. testing s&s=0.75',
    'base_case_lockdown_opt_with_limited_sens075_general_testing':  'LD, w. testing sens=0.75',
    'base_case_lockdown_opt_with_limited_spec075_general_testing':  'LD, w. testing spec=0.75',
    'base_case_lockdown_opt_with_limited_sens090_general_testing':  'LD, w. testing sens=0.90',
    'base_case_lockdown_opt_with_limited_spec090_general_testing':  'LD, w. testing spec=0.90',
    'base_case_lockdown_opt_with_limited_imperfect(0.85)_general_testing':  'LD, w. testing s&s=0.85',
    'base_case_lockdown_opt_with_limited_imperfect(0.90)_general_testing':  'LD, w. testing s&s=0.90',
    'base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing':  'LD, w. testing s&s=0.95',
    'base_case_lockdown_opt_3d_delay':  'LD, 3d delay',
    'base_case_lockdown_opt_7d_delay':  'LD, 7d delay',
    'base_case_lockdown_opt_14d_delay':  'LD, 14d delay',
    'base_case_lockdown_opt_28d_delay':  'LD, 28d delay',
    'base_case_6d_incubation': 'LD 6d incub.',
    'base_case_8d_incubation': 'LD 8d incub.',
    'test_and_trace_no_control_eta50': 'no control w CT, eta=50',
    'test_and_trace_no_control_eta100':'no control w CT, eta=100',
    'test_and_trace_no_control_R0_4.0_eta10': 'no control, R0=4.0, w CT, eta=10',
    'test_and_trace_no_control_R0_4.0_eta50': 'no control, R0=4.0, w CT, eta=50',
    'test_and_trace_no_control_R0_4.0_eta100': 'no control, R0=4.0, w CT, eta=100',
    'test_and_trace_lockdown_opt_eta10':    'LD w CT, eta=10, max 3M tests',
    'test_and_trace_lockdown_opt_eta50':    'LD w CT, eta=50, max 3M tests',
    'test_and_trace_lockdown_opt_eta75':    'LD w CT, eta=75, max 3M tests',
    'test_and_trace_lockdown_opt_eta95':    'LD w CT, eta=95, max 3M tests',
    'test_and_trace_lockdown_opt_eta100':    'LD w CT, eta=100, max 3M tests',
    'test_and_trace_lockdown_opt_eta50_R04':    'LD w CT, R0=4.0, eta=50, max 3M tests',
    'test_and_trace_lockdown_opt_eta75_R04':    'LD w CT, R0=4.0, eta=75, max 3M tests',
    'test_and_trace_lockdown_opt_eta100_R04':    'LD w CT, R0=4.0, eta=100, max 3M tests',
    'test_and_trace_lockdown_opt_eta50_R04_delta10':    'LD w CT, R0=4.0, eta=50, delta=10d, max 3M tests',
    'test_and_trace_lockdown_opt_eta75_R04_delta10':    'LD w CT, R0=4.0, eta=75, delta=10d, max 3M tests',
    'test_and_trace_lockdown_opt_eta100_R04_delta10':    'LD w CT, R0=4.0, eta=100, delta=10d, max 3M tests',
    'combo_base_case': 'combo, max 3M tests',
    'combo_base_case_R0_4.0': 'combo, R0=4.0, max 3M tests',
    'combo_sens_spec_0.95': 'combo, s&s=0.95, max 3M tests',
    'combo_sens_spec_0.85': 'combo, s&s=0.85, max 3M tests',
    'combo_R0_4.0_sens_spec_0.95': 'combo, R0=4.0, s&s=0.95, max 3M tests',
    'combo_R0_4.0_sens_spec_0.85': 'combo, R0=4.0, s&s=0.85, max 3M tests',
    'base_case_no_control_1000': 'no control, 1000 initial infected',
    'base_case_no_control_2000': 'no control, 2000 initial infected',
    'base_case_no_control_10000': 'no control, 10000 initial infected',
    'romer_tc_1000000_sens_spec_085': 'MT, s&s=0.85, max 1M tests',
    'romer_tc_2500000_sens_spec_085': 'MT, s&s=0.85, max 2.5M tests',
    'romer_tc_5000000_sens_spec_085': 'MT, s&s=0.85, max 5M tests',
    'romer_tc_25000000_sens_spec_085': 'MT, s&s=0.85, max 25M tests',
    'romer_tc_50000000_sens_spec_085': 'MT, s&s=0.85, max 50M tests',
    'romer_tc_100000000_sens_spec_085': 'MT, s&s=0.85, max 100M tests',
    'romer_tc_1000000': 'MT, max 1M tests',
    'romer_tc_2500000': 'MT, max 2.5M tests',
    'romer_tc_5000000': 'MT, max 5M tests',
    'romer_tc_25000000': 'MT, max 25M tests',
    'romer_tc_50000000': 'MT, max 50M tests',
    'romer_tc_100000000': 'MT, max 100M tests',
    'combo_base_case_test_and_trace': 'combo CT, max 3M tests',
    'combo_base_case_tc_1000000': 'combo max 1M tests',
    'combo_base_case_tc_2500000': 'combo max 2.5M tests',
    'combo_base_case_tc_5000000': 'combo max 5M tests',
    'combo_base_case_tc_25000000': 'combo max 25M tests',
    'combo_base_case_tc_50000000': 'combo max 50M tests',
    'combo_base_case_tc_100000000': 'combo max 100M tests',
    'combo_base_case_test_and_trace_tc1000000': 'combo max 1M tests + CT',
    'combo_base_case_test_and_trace_tc2500000': 'combo max 2.5M tests + CT',
    'combo_base_case_test_and_trace_tc5000000': 'combo max 5M tests + CT',
    'combo_base_case_test_and_trace_tc25000000': 'combo max 25M tests + CT',
    'combo_base_case_test_and_trace_tc50000000': 'combo max 50M tests + CT',
    'combo_base_case_test_and_trace_tc100000000': 'combo max 100M tests + CT',
    'combo_base_case_test_and_trace_ss085': 'combo max 3M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc1000000_ss085': 'combo max 1M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc2500000_ss085': 'combo max 2.5M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc5000000_ss085': 'combo max 5M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc25000000_ss085': 'combo max 25M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc50000000_ss085': 'combo max 50M tests + CT, s&s=0.85',
    'combo_base_case_test_and_trace_tc100000000_ss085': 'combo max 100M tests + CT, s&s=0.85',
    'combo_base_case_tc_50000000_sens_spec_085': 'combo max 50M tests, s&s=0.85',
    'combo_base_case_tc_25000000_sens_spec_085': 'combo max 25M tests, s&s=0.85',
    'combo_base_case_tc_5000000_sens_spec_085': 'combo max 5M tests, s&s=0.85',
    'combo_base_case_tc_2500000_sens_spec_085': 'combo max 2.50M tests, s&s=0.85',
    'combo_base_case_tc_1000000_sens_spec_085': 'combo max 1M tests, s&s=0.85',
    'combo_base_case_tc_100000000_sens_spec_085': 'combo max 100M tests, s&s=0.85',
    'combo_base_case_sens_spec_085': 'combo max 3M tests, s&s=0.85'

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
                axes[row,col].plot(time_steps, runs_data[run][outputs_array[row,col]], label=short_run_labels[run])
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

def pareto_plot(runs, obj_dict):

    #fig, axes = plt.subplots(ncols=2, figsize=(16, 8))
    fig, axes = plt.subplots(ncols=1, figsize=(8, 8))

    for run in runs:

        axes.scatter(obj_dict[run].iloc[:,0], -obj_dict[run].iloc[:,1], label=run_labels[run])
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


def sample_clouds(run, result_set, cols=2):
    medoid_solutions = {}
    medoid_obj = {}
    sample_obj = {}

    medoid_df = pd.read_csv('active_results/risk_analysis/' + run + '_' + result_set + '_risk.csv', delimiter=',')
    medoid_obj[run] = medoid_df[['Deaths', 'Economic impact']]
    medoid_solutions[run] = medoid_df.drop(columns=['Deaths', 'Economic impact'])

    n_sol = len(medoid_df.index)
    rows = int(np.ceil(n_sol / cols))

    color = iter(cm.tab10(np.linspace(0, 1, n_sol)))

    sample_obj[run] = {}

    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(20, 10))

    for row in range(0, rows):
        for col in range(0, cols):
            i = row * cols + col

            sample_df = pd.read_csv('active_results/risk_analysis/' + run + '__' + str(i) + '.csv', delimiter=',')
            sample_obj[run][i] = sample_df[['Deaths', 'Output']]

            index = medoid_df.index[i]

            c = next(color)
            axes[row, col].scatter(medoid_obj[run].Deaths[index], -medoid_obj[run]["Economic impact"][index], color=c)
            axes[row, col].scatter(sample_obj[run][i].Deaths, -sample_obj[run][i].Output, color=c, alpha=0.1, s=10)

    return fig