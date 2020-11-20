import numpy as np
from policy_epidemic_model_code import optimizable_corona_model

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
def create_epidemic_model(ksi_base=0,
               A_rel=0.5,
               r_AP=0,
               r_U=0.10,
               r_P=0.0,
               r_N=0.98,
               d_vaccine=800*14, # this is in _time steps_, not in days (days would be better)
               rel_rho=1.0,
               delta_param=5,
               omegaR_param=14,
               pii_D=0.01,
               R_0=2.5,
               rel_lambda_param=0.5,
               gamma_param=180.0,
               initial_infect=300,
               daily_testing_rate=0.0,
               testing_sensitivity=1.0,
               testing_specificity=1.0,
               tau_TT_daily=0.0,
               eta=0.0,
               unknown_q_rate=0.0,
               recovered_q_rate=0.0,
               negative_q_rate=0.0,
               positive_q_rate=0.999,
               testing_cost=100,
               ):

    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect, testing_cost, eta, gamma_param)

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

