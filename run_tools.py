import numpy as np
from policy_epidemic_model_code import optimizable_corona_model
from run_definitions import (ksi_base_default, A_rel_default, r_AP_default, r_U_default, r_P_default, r_N_default, d_vaccine_default, rel_rho_default, delta_param_default,
    omegaR_param_default, pii_D_default, R_0_default, rel_lambda_param_default, gamma_param_default, initial_infect_default, daily_testing_rate_default,
    testing_sensitivity_default, testing_specificity_default, tau_TT_daily_default, eta_default, unknown_q_rate_default, recovered_q_rate_default,
    negative_q_rate_default, positive_q_rate_default, testing_cost_default)

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
                    pop_size=60,
                    n_offsprings=30,
                    # sampling=get_sampling("real_random"),
                    # crossover=get_crossover("real_sbx", prob=0.9, eta=10),
                    # mutation=get_mutation("real_pm", eta=8),
                    # termination = get_termination("n_gen",50),

                    lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    lockdown_policy_lower_limits=list(0.5 * np.ones(15)),
                    # can't use len(l_p_c_d) within function param def
                    lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
                    testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
                    testing_policy_lower_limits=list(np.zeros(15)),
                    testing_policy_upper_limits=list(0.02 * np.ones(15)),
                    max_daily_tests=10000000,
                    p_ICU=0.01,
                    C_hos=100000,
                    T_rec=0.5  # recovery time in years from end of experiment
                    ):
    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                                     omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect, testing_cost, eta,
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

def create_optimization_run(
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
                   n_max_evals=100000,
               ),

               lockdown_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
               lockdown_policy_lower_limits=list(0.5 * np.ones(15)),  # can't use len(l_p_c_d) within function param def
               lockdown_policy_upper_limits=list(1.0 * np.ones(15)),  # needs to be different from lower limit
               testing_policy_control_days=[1, 15, 30, 60, 90, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600],
               testing_policy_lower_limits=list(np.zeros(15)),
               testing_policy_upper_limits=list(0.02 * np.ones(15)),
               max_daily_tests=10_000_000,
                p_ICU=0.01,
               C_hos=100000,
               T_rec=0.5, # recovery time in years from end of experiment
                **epidemic_model_params
               ):

    model, model_case = create_epidemic_model(epidemic_model_params)

    #print("DEBUG policy parameters:")
    #print("ld days:", lockdown_policy_control_days)
    #print("ld lolim: ", lockdown_policy_lower_limits)
    #print("ld hilim:", lockdown_policy_upper_limits)
    #print("t days:", testing_policy_control_days)
    #print("t lolim: ", testing_policy_lower_limits)
    #print("t hilim:", testing_policy_upper_limits)

    problem = COVID_policy(model, model_case, lockdown_policy_control_days, lockdown_policy_lower_limits,
                           lockdown_policy_upper_limits, testing_policy_control_days, testing_policy_lower_limits,
                           testing_policy_upper_limits, max_daily_tests, p_ICU, C_hos, T_rec)

    # create initial population here

    initial_pop_x = get_sampling("real_random") # random sampling for first population

    algorithm = NSGA2(
        pop_size=pop_size,
        n_offsprings=n_offsprings,
        sampling=initial_pop_x,
        crossover=crossover,
        mutation=mutation,
        eliminate_duplicates=True
    )

    return problem, algorithm, termination, model, model_case


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
