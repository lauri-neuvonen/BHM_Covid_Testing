# This script is used to run different optimization cases for different epidemic scenarios. It is a batch run compatible
# ...version of the notebook 'Covid_Policy_Optimization.ipynb'

# If run from command line, requires 2 arguments: number of max generations and list of runs to optimize for NSGA-II

# Dictionaries to hold run data and results:


import sys

import numpy as np
from pymoo.util.misc import stack
from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.optimize import minimize
# from pymoo.visualization.scatter import Scatter
from policy_epidemic_model_code import *
# import importlib
# import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume
from pymoo.util.termination.default import MultiObjectiveDefaultTermination
from pymoo.model.evaluator import Evaluator
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Run optimization run using an NSGA-II algorithm for COVID-19 simulator')
parser.add_argument('max_gen', type=int, help='maximum number of generations for NSGA-II algorithm.')
parser.add_argument('runs', type=str, nargs='+', help='dictionary of run values to be changed from defaults')

args = parser.parse_args()

epidemic_simulators = {}
result_dataframes = {}
objective_dataframes = {}
problems = {}

### RUN SETTINGS ###

max_gen = args.max_gen # 1000 # set low for testing, high for proper optimization runs. By default, optimization terminates
                    # ...when convergence has been observed or this limit of iterations reached.

# Tools for building optimization runs based on params.

arg_runs = args.runs # dictionary to hold all runs
print("Will optimize: ", arg_runs)

runs = {}

runs['base_case_no_control']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'lockdown_policy_control_days': [10000],   # no adjustments to testing policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'termination': get_termination("n_gen",1) # no optimization really...
}

runs['base_case_no_control_R0_4.0']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'lockdown_policy_control_days': [10000],   # no adjustments to testing policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'R_0': 4.0,
    'termination': get_termination("n_gen",1) # no optimization really...
}
### ROMER CASE SCENARIOS ###
#------------------------------------------#

runs['romer']={
    'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0],
    'lockdown_policy_upper_limits': [0]
}

#------------------------------------------#

runs['romer_R0_4.0']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'R_0': 4.0, # set R0 to a higher value
}

runs['romer_R0_1.25']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'R_0': 1.25, # set R0 to a higher value
}

# #------------------------------------------#

runs['romer_R0_4.0_sens_spec_075']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.75,
    'testing_specificity': 0.75,
    'R_0': 4.0, # set R0 to a higher value
}


runs['romer_R0_4.0_sens_spec_085']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.85,
    'testing_specificity': 0.85,
    'R_0': 4.0, # set R0 to a higher value
}

runs['romer_R0_4.0_sens_spec_090']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.90,
    'testing_specificity': 0.90,
    'R_0': 4.0, # set R0 to a higher value
}

runs['romer_R0_4.0_sens_spec_095']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.95,
    'testing_specificity': 0.95,
    'R_0': 4.0, # set R0 to a higher value
}
#------------------------------------------#

runs['romer_6d_incubation']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'delta_param': 6
}

#------------------------------------------#

runs['romer_8d_incubation']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'delta_param': 8
}

#------------------------------------------#

runs['romer_8d_incubation_sens_spec_075']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'delta_param': 8,
    'testing_sensitivity': 0.75,
    'testing_specificity': 0.75
}

#------------------------------------------#

runs['romer_8d_incubation_sens_spec_090']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'delta_param': 8,
    'testing_sensitivity': 0.75,
    'testing_specificity': 0.75
}


runs['romer_3d_delay']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_policy_control_days': [3, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
}

runs['romer_7d_delay']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_policy_control_days': [7, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
}

runs['romer_14d_delay']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_policy_control_days': [14, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
}

runs['romer_28d_delay']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_policy_control_days': [28, 29, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
}

runs['romer_sens_075']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.75,
}

runs['romer_spec_075']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_specificity': 0.75,
}

runs['romer_sens_spec_075']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.75,
    'testing_specificity': 0.75,
}

#------------------------------------------#

runs['romer_sens_spec_085']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.85,
    'testing_specificity': 0.85,
}

#------------------------------------------#

runs['romer_sens_spec_090']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.90,
    'testing_specificity': 0.90,
}

#------------------------------------------#

runs['romer_sens_spec_095']={
    'lockdown_policy_control_days': [10000],   # no adjustments to lockdown policy
    'lockdown_policy_lower_limits': [0.0],
    'lockdown_policy_upper_limits': [0.05],
    'testing_sensitivity': 0.95,
    'testing_specificity': 0.95,
}

#------------------------------------------#


### LOCKDOWN OPTIMIZATION CASE SCENARIOS ###

#------------------------------------------#

runs['base_case_lockdown_opt']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

#------------------------------------------#

runs['base_case_lockdown_opt_R0_4.0']={
    'R_0': 4.0,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

#------------------------------------------#

runs['base_case_lockdown_opt_with_limited_general_testing']={
    'testing_rate': 0.005,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

#------------------------------------------#

runs['base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.75,
    'testing_specificity': 0.75,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_with_limited_sens075_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.75,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_with_limited_spec075_general_testing']={
    'testing_rate': 0.005,
    'testing_specificity': 0.75,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_with_limited_sens090_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.90,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_with_limited_spec090_general_testing']={
    'testing_rate': 0.005,
    'testing_specificity': 0.90,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}



runs['base_case_lockdown_opt_with_limited_imperfect(0.85)_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.85,
    'testing_specificity': 0.85,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}
runs['base_case_lockdown_opt_with_limited_imperfect(0.90)_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.90,
    'testing_specificity': 0.90,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}
runs['base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing']={
    'testing_rate': 0.005,
    'testing_sensitivity': 0.95,
    'testing_specificity': 0.95,
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}
#------------------------------------------#

runs['base_case_lockdown_opt_3d_delay']={
    'lockdown_policy_control_days': [3, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_7d_delay']={
    'lockdown_policy_control_days': [7, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_14d_delay']={
    'lockdown_policy_control_days': [14, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_lockdown_opt_28d_delay']={
    'lockdown_policy_control_days': [28, 29, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05]
}

runs['base_case_6d_incubation']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'delta_param': 6,
}


runs['base_case_8d_incubation']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'delta_param': 8,
}

##### TEST and TRACE CASES #####


# Fixed testing rate:
runs['test_and_trace_lockdown_opt_eta50']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.50,
    'tau_TT': 0.2,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta75']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.75,
    'tau_TT': 0.2,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta95']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.95,
    'tau_TT': 0.2,
    'r_U': 0.01
}


runs['test_and_trace_lockdown_opt_eta100']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 1.00,
    'tau_TT': 0.2,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta50_R04']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.50,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'r_U': 0.01
}
runs['test_and_trace_lockdown_opt_eta75_R04']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.75,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta100_R04']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 1.00,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta50_R04_delta10']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.50,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'delta_param': 10,
    'r_U': 0.01
}

runs['test_and_trace_lockdown_opt_eta75_R04_delta10']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 0.75,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'delta_param': 10,
    'r_U': 0.01
}


runs['test_and_trace_lockdown_opt_eta100_R04_delta10']={
    'testing_policy_control_days': [10000],   # no adjustments to testing policy
    'testing_policy_lower_limits': [0.0],
    'testing_policy_upper_limits': [0.05],
    'eta': 1.00,
    'R_0': 4.0,
    'tau_TT': 0.2,
    'delta_param': 10,
    'r_U': 0.01
}


#### OPTIMIZATION ENGINE ###
# MODIFY ONLY IF YOU KNOW WHAT YOU'RE _DOING
# Modifications here affect how the optimizer works, not run definitions

class COVID_policy(Problem):

    def __init__(self, model, model_case, lockdown_policy_control_days, lockdown_policy_lower_limits,
                 lockdown_policy_upper_limits, testing_policy_control_days, testing_policy_lower_limits,
                 testing_policy_upper_limits):
        self.model = model
        self.model_case = model_case


        if lockdown_policy_control_days == "NA":
            n_var_ld = 0
            self.lockdown_policy_control_days = []
            self.lockdown_policy_lower_limits = []
            self.lockdown_policy_upper_limits = []
        else:
            n_var_ld = len(lockdown_policy_control_days)
            self.lockdown_policy_control_days = lockdown_policy_control_days
            self.lockdown_policy_upper_limits = lockdown_policy_upper_limits
            self.lockdown_policy_lower_limits = lockdown_policy_lower_limits

        if testing_policy_control_days == "NA":
            n_var_t = 0
            self.testing_policy_lower_limits = []
            self.testing_policy_upper_limits = []
        else:
            n_var_t = len(testing_policy_control_days)
            self.testing_policy_control_days = testing_policy_control_days
            self.testing_policy_lower_limits = testing_policy_lower_limits
            self.testing_policy_upper_limits = testing_policy_upper_limits


        super().__init__(n_var=n_var_ld+n_var_t,
                         n_obj=2,
                         n_constr=0,
                         xl=np.array(self.lockdown_policy_lower_limits + self.testing_policy_lower_limits),
                         xu=np.array(self.lockdown_policy_upper_limits + self.testing_policy_upper_limits)
        )

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = []
        f2 = []
        f3 = []
        for j in range(len(x[:, 1])):  # evaluate f1 and f2 for each individual

            # these slices contain indexing for decision variables representing lockdown and testing:
            lockdown_var_slice = slice(0, len(self.lockdown_policy_control_days) + 1)
            testing_var_slice = slice(len(self.lockdown_policy_control_days),
                                      len(self.lockdown_policy_control_days) + len(
                                          self.testing_policy_control_days) + 1)

            # create policies (dictionaries) for lockdown and testing
            lockdown_policy = create_sub_policy(self.lockdown_policy_control_days, x[j, lockdown_var_slice])
            testing_policy = create_sub_policy(self.testing_policy_control_days, x[j, testing_var_slice])
            policy = Policy(lockdown_policy, testing_policy)

            Reported_D, Notinfected_D, Unreported_D, Infected_D, \
            False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total, total_cost, tests, Unk_NA_nQ_D, Unk_NA_Q_D, K_NA_nQ_D, Unk_IA_nQ_D, Unk_IA_Q_D, K_IA_Q_D, alpha_D, ksi_TT_D, Symptomatic_D \
                = self.model.solve_case(self.model_case, policy)

            # objectives scaled to roughly same scale
            f1.append(Dead_D[-1] * self.model.pop / 1000)
            f2.append(-Y_total / (14 * 365 * self.model.T_years))
            f3.append(total_cost / (self.model.pop * Y_total / (
                        14 * 365 * self.model.T_years)))  # algorithm minimizes total cost per (saved) output unit

            # f1.append(Dead_D[-1])
            # f2.append(-Y_total) # algorithm minimizes, Y_total needs to be max'd -> negative

        out["F"] = np.column_stack([f1, f2, f3])
        # out["G"] = np.column_stack([g1, g2])


class Policy():
    def __init__(self, lockdown_policy, testing_policy):
        self.lockdown_policy = lockdown_policy
        self.testing_policy = testing_policy


def create_sub_policy(policy_control_times, policy_control_values):
    policy = {}  # this will hold the policy in format suitable for input to the epidemic model
    if policy_control_times == "NA":
        return "NA"
    else:
        for (i, t) in enumerate(policy_control_times):
            policy[t] = policy_control_values[i]

        return policy

def create_policy(lockdown_policy, testing_policy):

    return Policy(lockdown_policy, testing_policy)

# Run generator

# NOTE: default values for all adjustable run parameters defined in function definition below:
def create_run(ksi_base=0,
               A_rel=0.5,
               r_AP=0,
               r_U=0.10,
               d_vaccine=800,
               rel_rho=1.0,
               delta_param=5,
               omegaR_param=14,
               pii_D=0.01,
               R_0=2.5,
               rel_lambda_param=0.5,
               gamma_param=180.0,
               initial_infect=300,
               testing_rate=0.0,
               testing_sensitivity=1.0,
               testing_specificity=1.0,
               tau_TT=0.5,
               eta=0.0,
               unknown_q_rate=0.0,
               recovered_q_rate=0.0,
               negative_q_rate=0.0,
               positive_q_rate=0.999,
               testing_cost=100,
               pop_size=60,
               n_offsprings=30,
               sampling=get_sampling("real_random"),
               crossover=get_crossover("real_sbx", prob=0.9, eta=10),
               mutation=get_mutation("real_pm", eta=8),
               eliminate_duplicates=True,
               filename="foo",
               # termination = get_termination("n_gen",50),
               termination=MultiObjectiveDefaultTermination(
                   x_tol=1e-8,
                   cv_tol=1e-6,
                   f_tol=0.0025,
                   nth_gen=5,
                   n_last=30,
                   n_max_gen=max_gen,
                   n_max_evals=100000
               ),

               lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
               lockdown_policy_lower_limits=list(0.5 * np.ones(15)),  # can't use len(l_p_c_d) within function param def
               lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
               testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
               testing_policy_lower_limits=list(np.zeros(15)),
               testing_policy_upper_limits=list(0.2 * np.ones(15))
               ):
    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect, testing_cost, eta, gamma_param)

    model_case = {
        'tau_paramA': testing_rate,
        'test_sens': testing_sensitivity,
        'test_spec': testing_specificity,
        'tau_TT': tau_TT,
        'ksi_U': (1 + unknown_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_P': (1 + positive_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_N': (1 + negative_q_rate) ** (1. / model.Delta_time) - 1,
        'ksi_R': (1 + recovered_q_rate) ** (1. / model.Delta_time) - 1,
        'r_U': r_U,  # should be redundant!
        'r_P': (1 + 0.98) ** (1. / model.Delta_time) - 1,
        'r_AP': 0.0,
        'r_N': (1 + 0.98) ** (1. / model.Delta_time) - 1,
        'r_R': (1 + 0.999) ** (1. / model.Delta_time) - 1,
        'd_start_exp': 0.,
        'experiment': "baseline_vaccine_tag"
    }
    print("DEBUG policy parameters:")
    print("ld days:", lockdown_policy_control_days)
    print("ld lolim: ", lockdown_policy_lower_limits)
    print("ld hilim:", lockdown_policy_upper_limits)
    print("t days:", testing_policy_control_days)
    print("t lolim: ", testing_policy_lower_limits)
    print("t hilim:", testing_policy_upper_limits)

    problem = COVID_policy(model, model_case, lockdown_policy_control_days, lockdown_policy_lower_limits,
                           lockdown_policy_upper_limits, testing_policy_control_days, testing_policy_lower_limits,
                           testing_policy_upper_limits)

    # create initial population here
    try:
        if (len(lockdown_policy_control_days) > 1) and (len(testing_policy_control_days) == 1):
            initial_pop_x = pd.read_csv('results/base_case_lockdown_opt_results.csv', delimiter=',').to_numpy()
            initial_pop = Evaluator().eval(problem, initial_pop_x)
        elif (len(lockdown_policy_control_days) == 1) and (len(testing_policy_control_days) < 1):
            initial_pop_x = pd.read_csv('romer_results.csv', delimiter=',').to_numpy()
            initial_pop = Evaluator().eval(problem, initial_pop_x)
        else:
            initial_pop_x = get_sampling("real_random")
    except:
        initial_pop_x = get_sampling("real_random")

    algorithm = NSGA2(
        pop_size=pop_size,
        n_offsprings=n_offsprings,
        sampling=initial_pop_x,
        crossover=crossover,
        mutation=mutation,
        eliminate_duplicates=True
    )

    return problem, algorithm, termination, model, model_case


# loop through the different runs and
for run in arg_runs:
    problem, algorithm, termination, model, model_case = create_run(**runs[run])
    epidemic_simulators[run] = (model, model_case)
    problems[run] = problem

    print("Optimizing run ", run)
    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=1,
                   # pf=problem.pareto_front(use_cache=False),
                   save_history=False,  # True Needed for convergence analysis
                   verbose=True)

    # Which one to use? Pandas makes it easier to save column names...
    # np.savetxt('results/'+run+'_results.csv', res.X, delimiter=",")
    # np.savetxt('results/'+run+'_objectives.csv', res.F, delimiter=",")

    res_df = pd.DataFrame(data=res.X,
                          columns=problem.lockdown_policy_control_days + problem.testing_policy_control_days)
    res_df.to_csv('results/' + run + '_results.csv', index=False)
    result_dataframes[run] = res_df

    obj_df = pd.DataFrame(data=res.F, columns=['Deaths', 'Output', 'Direct Cost'])
    obj_df.to_csv('results/' + run + '_objectives.csv', index=False)
    objective_dataframes[run] = obj_df

print("Runs completed and results saved")

