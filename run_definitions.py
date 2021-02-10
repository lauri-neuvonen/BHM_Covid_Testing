
from pymoo.factory import get_termination
import numpy as np
from pymoo.factory import get_sampling, get_crossover, get_mutation


### DEFAULT VALUES ###

# NOTE! these names should be unique to allow easier import with from ... import *

### Epidemic model ###

ksi_base_default = 0
A_rel_default = 0.5
r_AP_default = 0
r_U_default = 0.10
r_P_default = 0.0
r_N_default = 0.98
d_vaccine_default = 800*14 # this is in _time steps_, not in days (days would be better)
rel_rho_default = 1.0
delta_param_default = 6.71
omegaR_param_default = 14
pii_D_default = 0.01
R_0_default = 2.5
rel_lambda_param_default = 0.5
lambda_param_default = 1.0 # NOTE: Currenly a normalized constant. Changing this might have unexpected consequences.
gamma_param_default = 180.0
initial_infect_default = 20000 # testing number based on no control tests. 0.2 per mille of total population
daily_testing_rate_default = 0.0
testing_sensitivity_default = 1.0
testing_specificity_default = 1.0
tau_TT_daily_default = 0.0
eta_default = 0.0
unknown_q_rate_default = 0.0
recovered_q_rate_default = 0.0
negative_q_rate_default = 0.0
positive_q_rate_default = 0.999
testing_cost_default = 100

### Optimization ###

# policy setup

lockdown_policy_control_days_def = [1, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570]
lockdown_policy_lower_limits_def = list(0.5 * np.ones(20))
lockdown_policy_upper_limits_def = list(1.0 * np.ones(20))
testing_policy_control_days_def = [1, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570]
testing_policy_lower_limits_def = list(np.zeros(20))
testing_policy_upper_limits_def = list(0.04 * np.ones(20))
max_daily_tests_def = 3000000
p_ICU_def = 0.01
C_hos_def = 30000 # ICU capacity
T_rec_def = 0.5

# NSGA-II parameters

pop_size_def=60
n_offsprings_def=30
sampling_def=get_sampling("real_random")
crossover_def=get_crossover("real_sbx", prob=0.9, eta=10)
mutation_def=get_mutation("real_pm", eta=8)

x_tol_def=1e-8
cv_tol_def=1e-6
f_tol_def=0.0025
nth_gen_def=5
n_last_def=30
n_max_evals_def=100000

def get_runs_definitions():

    # runs are defined below as modifications to the default parameter values. If parameter name is not mentioned
    # in the run definition, its default value will be used. Default values are defined in Covid_Run_Optimizer.py
    # (and possibly separately elsewhere, like in visualizations).

    runs = {}

    runs['base_case_no_control'] = {
        'testing_policy_control_days': "NA",  # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'lockdown_policy_control_days': "NA",  # no adjustments to testing policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'termination': get_termination("n_gen", 1)  # no optimization really...
    }

    runs['base_case_no_control_R0_4.0']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'lockdown_policy_control_days': "NA",   # no adjustments to testing policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'R_0': 4.0,
        'termination': get_termination("n_gen",1) # no optimization really...
    }
    ### ROMER CASE SCENARIOS ###
    #------------------------------------------#

    runs['romer']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': []
    }

    runs['romer_ICUC_50000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'C_hos': 50000
    }

    runs['romer_ICUC_10000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'C_hos': 10000
    }

    runs['romer_ICUC_5000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'C_hos': 5000
    }

    runs['romer_ICUC_1000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'C_hos': 1000
    }

    runs['romer_no_limit']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_policy_upper_limits': list(0.05 * np.ones(15)),
        'max_daily_tests': 100000000,
    }

    runs['romer_terminal_0.25']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'T_rec': 0.25
    }

    runs['romer_terminal_1.0']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'T_rec': 1.0
    }

    #------------------------------------------#

    runs['romer_R0_4.0']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'R_0': 4.0, # set R0 to a higher value
    }

    runs['romer_R0_4.0_no_limit']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'R_0': 4.0, # set R0 to a higher value
        'max_daily_tests': 100000000,
    }

    runs['romer_R0_1.25']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'R_0': 1.25, # set R0 to a lower value
    }

    # #------------------------------------------#

    runs['romer_R0_4.0_sens_spec_075']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.75,
        'testing_specificity': 0.75,
        'R_0': 4.0, # set R0 to a higher value
    }


    runs['romer_R0_4.0_sens_spec_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'R_0': 4.0, # set R0 to a higher value
    }

    runs['romer_R0_4.0_sens_spec_090']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.90,
        'testing_specificity': 0.90,
        'R_0': 4.0, # set R0 to a higher value
    }

    runs['romer_R0_4.0_sens_spec_095']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.95,
        'testing_specificity': 0.95,
        'R_0': 4.0, # set R0 to a higher value
    }
    #------------------------------------------#

    runs['romer_6d_incubation']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'delta_param': 6
    }

    runs['romer_6d_incubation_sens_spec_090']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'delta_param': 6,
        'testing_sensitivity': 0.90,
        'testing_specificity': 0.90
    }

    #------------------------------------------#

    runs['romer_8d_incubation']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'delta_param': 8
    }

    #------------------------------------------#

    runs['romer_8d_incubation_sens_spec_075']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'delta_param': 8,
        'testing_sensitivity': 0.75,
        'testing_specificity': 0.75
    }

    #------------------------------------------#

    runs['romer_8d_incubation_sens_spec_090']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'delta_param': 8,
        'testing_sensitivity': 0.90,
        'testing_specificity': 0.90
    }


    runs['romer_3d_delay']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_policy_control_days': [3, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
    }


    runs['romer_7d_delay']={    # NOT VERY INFORMATIVE...
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_policy_control_days': [7, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
    }


    runs['romer_14d_delay']={   # NOT VERY INFORMATIVE...
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_policy_control_days': [14, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
    }

    runs['romer_28d_delay']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_policy_control_days': [28, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
    }

    runs['romer_sens_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
    }

    runs['romer_spec_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_specificity': 0.85,
    }


    #------------------------------------------#

    runs['romer_sens_spec_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['romer_sens_spec_085_low_mut_param'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'mutation': get_mutation("real_pm", eta=4)
    }

    runs['romer_sens_spec_085_hi_mut_param'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'mutation': get_mutation("real_pm", eta=12)
    }

    runs['romer_sens_spec_085_low_cross_param'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'crossover': get_crossover("real_sbx", prob=0.9, eta=5)
    }

    runs['romer_sens_spec_085_hi_cross_param'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'crossover': get_crossover("real_sbx", prob=0.9, eta=15)
    }

    #------------------------------------------#

    runs['romer_sens_spec_090']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.90,
        'testing_specificity': 0.90,
    }

    #------------------------------------------#

    runs['romer_sens_spec_095']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.95,
        'testing_specificity': 0.95,
    }

    runs['romer_sens_spec_098']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.98,
        'testing_specificity': 0.98,
    }

    runs['romer_sens_spec_099']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.99,
        'testing_specificity': 0.99,
    }

    runs['romer_tc_500000_sens_spec_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 500000
    }

    runs['romer_tc_1000000_sens_spec_085']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 1000000
    }

    runs['romer_tc_2500000_sens_spec_085'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 2500000
    }

    runs['romer_tc_5000000_sens_spec_085'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 5000000
    }

    runs['romer_tc_25000000_sens_spec_085'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 25000000
    }

    runs['romer_tc_50000000_sens_spec_085'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 50000000
    }

    runs['romer_tc_100000000_sens_spec_085'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 100000000
    }

    runs['romer_tc_500000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 500000
    }

    runs['romer_tc_1000000']={
        'lockdown_policy_control_days': "NA",   # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 1000000
    }

    runs['romer_tc_2500000'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 2500000
    }

    runs['romer_tc_5000000'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 5000000
    }

    runs['romer_tc_25000000'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 25000000
    }

    runs['romer_tc_50000000'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 50000000
    }

    runs['romer_tc_100000000'] = {
        'lockdown_policy_control_days': "NA",  # no adjustments to lockdown policy
        'lockdown_policy_lower_limits': [],
        'lockdown_policy_upper_limits': [],
        'max_daily_tests': 100000000
    }

    #------------------------------------------#


    ### LOCKDOWN OPTIMIZATION CASE SCENARIOS ###

    #------------------------------------------#

    runs['base_case_lockdown_opt']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_ICUC_1000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'C_hos': 1000
    }

    runs['base_case_lockdown_opt_ICUC_5000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'C_hos': 5000
    }

    runs['base_case_lockdown_opt_ICUC_10000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'C_hos': 10000
    }

    runs['base_case_lockdown_opt_ICUC_50000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'C_hos': 50000
    }

    runs['base_case_lockdown_opt_terminal_0.25']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'T_rec': 0.25
    }

    runs['base_case_lockdown_opt_terminal_1.0']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'T_rec': 1
    }

    #------------------------------------------#

    runs['base_case_lockdown_opt_R0_4.0']={
        'R_0': 4.0,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    #------------------------------------------#

    runs['base_case_lockdown_opt_with_limited_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    #------------------------------------------#

    runs['base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.75,
        'testing_specificity': 0.75,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_with_limited_sens075_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.75,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_with_limited_spec075_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_specificity': 0.75,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_with_limited_sens090_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.90,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_with_limited_spec090_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_specificity': 0.90,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }



    runs['base_case_lockdown_opt_with_limited_imperfect(0.85)_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }
    runs['base_case_lockdown_opt_with_limited_imperfect(0.90)_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.90,
        'testing_specificity': 0.90,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }
    runs['base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing']={
        'daily_testing_rate': 0.01,
        'testing_sensitivity': 0.95,
        'testing_specificity': 0.95,
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }
    #------------------------------------------#

    runs['base_case_lockdown_opt_3d_delay']={
        'lockdown_policy_control_days': [3, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_7d_delay']={
        'lockdown_policy_control_days': [7, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_14d_delay']={
        'lockdown_policy_control_days': [14, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_lockdown_opt_28d_delay']={
        'lockdown_policy_control_days': [28, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570],
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': []
    }

    runs['base_case_6d_incubation']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'delta_param': 6,
    }


    runs['base_case_8d_incubation']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'delta_param': 8,
    }

    ##### TEST and TRACE CASES #####


    # Fixed testing rate:
    runs['test_and_trace_lockdown_opt_eta10']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.50,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }

    runs['test_and_trace_lockdown_opt_eta50']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.50,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }

    runs['test_and_trace_lockdown_opt_eta75']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.75,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000

    }

    runs['test_and_trace_lockdown_opt_eta95']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.95,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }


    runs['test_and_trace_lockdown_opt_eta100']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }

    runs['test_and_trace_lockdown_opt_eta100_ICUC_1000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'C_hos': 1000
    }

    runs['test_and_trace_lockdown_opt_eta100_ICUC_5000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'C_hos': 5000
    }

    runs['test_and_trace_lockdown_opt_eta100_ICUC_10000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'C_hos': 10000
    }

    runs['test_and_trace_lockdown_opt_eta100_ICUC_50000']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'C_hos': 50000
    }

    runs['test_and_trace_lockdown_opt_eta50_R04']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.50,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }
    runs['test_and_trace_lockdown_opt_eta75_R04']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.75,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }

    runs['test_and_trace_lockdown_opt_eta100_R04']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000
    }

    runs['test_and_trace_lockdown_opt_eta50_R04_delta10']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.50,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'delta_param': 10
    }

    runs['test_and_trace_lockdown_opt_eta75_R04_delta10']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 0.75,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'delta_param': 10
    }


    runs['test_and_trace_lockdown_opt_eta100_R04_delta10']={
        'testing_policy_control_days': "NA",   # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
        'eta': 1.00,
        'R_0': 4.0,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'delta_param': 10
    }

     #------- COMBO STRATEGIES ------- #

    runs['combo_base_case']={
    }

    runs['combo_base_case_ICUC_50000']={
        'C_hos':50000
    }

    runs['combo_base_case_ICUC_10000']={
        'C_hos':10000
    }

    runs['combo_base_case_ICUC_5000']={
        'C_hos':5000
    }

    runs['combo_base_case_ICUC_1000']={
        'C_hos':1000
    }

    runs['combo_base_case_tc_500000'] = {
        'max_daily_tests': 500000
    }
    runs['combo_base_case_tc_1000000'] = {
        'max_daily_tests': 1000000
    }

    runs['combo_base_case_tc_2500000'] = {
        'max_daily_tests': 2500000
    }

    runs['combo_base_case_tc_5000000'] = {
        'max_daily_tests': 5000000
    }

    runs['combo_base_case_tc_25000000'] = {
        'max_daily_tests': 25000000
    }

    runs['combo_base_case_tc_50000000'] = {
        'max_daily_tests': 50000000
    }

    runs['combo_base_case_tc_100000000'] = {
        'max_daily_tests': 100000000
    }

    runs['combo_base_case_sens_spec_085'] = {
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_base_case_tc_500000_sens_spec_085'] = {
        'max_daily_tests': 500000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }
    runs['combo_base_case_tc_1000000_sens_spec_085'] = {
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
        'max_daily_tests': 1000000
    }

    runs['combo_base_case_tc_2500000_sens_spec_085'] = {
        'max_daily_tests': 2500000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85
    }

    runs['combo_base_case_tc_5000000_sens_spec_085'] = {
        'max_daily_tests': 5000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85
    }

    runs['combo_base_case_tc_25000000_sens_spec_085'] = {
        'max_daily_tests': 25000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85
    }

    runs['combo_base_case_tc_50000000_sens_spec_085'] = {
        'max_daily_tests': 50000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85
    }

    runs['combo_base_case_R0_4.0'] = {
        'R_0': 4.0,
    }

    runs['combo_base_case_test_and_trace'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 10000000,
    }
    runs['combo_base_case_test_and_trace_tc500000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 500000,
    }
    runs['combo_base_case_test_and_trace_tc1000000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 1000000,
    }
    runs['combo_base_case_test_and_trace_tc2500000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 2500000,
    }

    runs['combo_base_case_test_and_trace_tc5000000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 5000000,
    }


    runs['combo_base_case_test_and_trace_tc25000000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 25000000,
    }

    runs['combo_base_case_test_and_trace_tc50000000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 50000000,
    }
    runs['combo_base_case_test_and_trace_tc100000000'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
    }

    runs['combo_base_case_test_and_trace_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_base_case_test_and_trace_tc100000000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 100000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }
    runs['combo_base_case_test_and_trace_tc50000000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 50000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_base_case_test_and_trace_tc25000000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 25000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_base_case_test_and_trace_tc5000000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 5000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_base_case_test_and_trace_tc1000000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 1000000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }


    runs['combo_base_case_test_and_trace_tc2500000_ss085'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'max_daily_tests': 2500000,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }





    runs['combo_lockdown_test_and_trace'] = {
        'eta': 1.00,
        'tau_TT_daily': 0.5,
        'r_U': 0.01,
        'testing_policy_control_days': "NA",  # no adjustments to testing policy
        'testing_policy_lower_limits': [],
        'testing_policy_upper_limits': [],
    }

    runs['combo_sens_spec_0.95'] = {
        'testing_sensitivity': 0.95,
        'testing_specificity': 0.95,
    }

    runs['combo_sens_spec_0.85'] = {
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    runs['combo_R0_4.0_sens_spec_0.95'] = {
        'R_0': 4.0,
        'testing_sensitivity': 0.95,
        'testing_specificity': 0.95,
    }

    runs['combo_R0_4.0_sens_spec_0.85'] = {
        'R_0': 4.0,
        'testing_sensitivity': 0.85,
        'testing_specificity': 0.85,
    }

    return runs
