import numpy as np
from multiprocessing import Pool
from itertools import product
import os

try:
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
    from plotly.offline import init_notebook_mode, iplot
except ImportError:
    print("Installing plotly. This may take a while.")
    from pip._internal import main as pipmain
    pipmain(['install', 'plotly'])
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
    from plotly.offline import init_notebook_mode, iplot

class optimizable_corona_model(object):
    # ksi_base: baseline quarantine rate
    # A_rel: relative productivity of quarantined
    # d_vaccine: date of vaccine
    # rel_rho: relative infectiousness?
    # delta_param: 1/ mean days to show symptoms
    # omegaR_param ~mean days to recovery (check details)
    # pii_D: mortality rate
    # R_0: basic reproduction number
    # rel_lambda_param: effectiveness of quarantine?
    # initial_infect

    def __init__(self, ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                 omegaR_param, pii_D, R_0, rel_lambda_paramQ, initial_infect, test_cost, eta, gamma_param):
        self.pop        = 340_000_000
        self.T_years    = 2
        self.Delta_time     = 14
        self.T          = self.T_years * 365 * self.Delta_time
        self.lambda_param          = 1
        self.gamma          =   1/(self.Delta_time*gamma_param)       # immunity loss rate
        self.ksi_base_high= .999
        self.r_high     = .999
        self.r          = .98
        self.r_AP       = r_AP
        self.delta          = 1/(self.Delta_time*delta_param)
        self.omegaR         = 1/(self.Delta_time*omegaR_param)

        self.omegaD         = self.omegaR*pii_D/(1-pii_D)
        self.lambda_paramQ         = rel_lambda_paramQ*self.lambda_param             # lambda for quarantine
        self.rhoS         = R_0/((self.lambda_param/self.delta)*(rel_rho + self.delta/(self.omegaR+self.omegaD)))
        self.rhoA         = rel_rho*self.rhoS

        self.InitialInfect = initial_infect
        self.d_vaccine     = d_vaccine
        self.A_rel         = A_rel
        self.ksi_base        = ksi_base
        self.test_cost     = test_cost

        self.sigma = 1/7
        self.sigma_Q = 0.00

        self.eta = eta # test and trace efficiency parameter

        #self.test_sens      = test_sens
        #self.test_spec      = test_spec

        self.baseline = {
            'tau_paramA'            : 0.,
            'test_sens'     : 1.0,
            'test_spec'     : 1.0,
            'ksi_U'           : 0., # baseline quarantine rate
            'ksi_P'           : 0.,
            'ksi_N'           : 0.,
            'ksi_R'           : 0.,
            'r_U'           : self.r,
            'r_P'           : 0.,
            'r_AP'          : self.r_AP,
            'r_N'           : self.r,

            'r_R'           : self.r_high,
            'd_start_exp'   : 0.,
            'experiment'    : "baseline_vaccine_tag"
        }

        tau_param_A_daily_target = 0
        ksi_U_daily_target = ksi_base
        ksi_P_daily_target = self.ksi_base_high
        ksi_N_daily_target = ksi_base
        ksi_R_daily_target = 0

        r_U_daily_target = 0
        r_N_daily_target = 0
        r_P_daily_target = 0
        r_AP_daily_target = self.r_AP  # self.r
        r_R_daily_target = self.r_high

        self.policy_offset = 14

        self.optimization_no_test = {
            'tau_paramA'            : 0,
            'test_sens'     : 1.0,
            'test_spec'     : 1.0,
            'ksi_U'           : 0,
            'ksi_P'           : (1+ksi_P_daily_target)**(1./self.Delta_time)-1,
            'ksi_N'           : 0,
            'ksi_R'           : 0,
            'r_U'           : (1+r_U_daily_target)**(1./self.Delta_time)-1, # should be redundant!
            'r_P'           : 0,
            'r_AP'          : 0,
            'r_N'           : (1+r_N_daily_target)**(1./self.Delta_time)-1,
            'r_R'           : (1+r_R_daily_target)**(1./self.Delta_time)-1,
            'd_start_exp': 0.,
            'experiment': "baseline_vaccine_tag"
        }

        self.optimization_general_testing_no_lockdown = {
            'tau_paramA': (1 + tau_param_A_daily_target) ** (1. / self.Delta_time) - 1,
            'test_sens': 1.0,
            'test_spec': 1.0,
            'ksi_U': 0,
            'ksi_P': (1 + ksi_P_daily_target) ** (1. / self.Delta_time) - 1,
            'ksi_N': 0,
            'ksi_R': 0,
            'r_U': (1 + r_U_daily_target) ** (1. / self.Delta_time) - 1,
            'r_P': 0,
            'r_AP': 0,
            'r_N': (1 + r_N_daily_target) ** (1. / self.Delta_time) - 1,
            'r_R': (1 + r_R_daily_target) ** (1. / self.Delta_time) - 1,
            'd_start_exp': 0.,
            'experiment': "baseline_vaccine_tag"
        }

        self.optimization_general_imperf_testing_no_lockdown = {
            'tau_paramA': (1 + tau_param_A_daily_target) ** (1. / self.Delta_time) - 1,
            'test_sens': 0.75,
            'test_spec': 0.75,
            'ksi_U': 0,
            'ksi_P': (1 + ksi_P_daily_target) ** (1. / self.Delta_time) - 1,
            'ksi_N': 0,
            'ksi_R': 0,
            'r_U': (1 + r_U_daily_target) ** (1. / self.Delta_time) - 1,
            'r_P': 0,
            'r_AP': 0,
            'r_N': (1 + r_N_daily_target) ** (1. / self.Delta_time) - 1,
            'r_R': (1 + r_R_daily_target) ** (1. / self.Delta_time) - 1,
            'd_start_exp': 0.,
            'experiment': "baseline_vaccine_tag"
        }

        self.optimization_general_testing_with_lockdown = {
            'tau_paramA': (1 + tau_param_A_daily_target) ** (1. / self.Delta_time) - 1,
            'test_sens': 1.0,
            'test_spec': 1.0,
            'ksi_U': 0,
            'ksi_P': (1 + ksi_P_daily_target) ** (1. / self.Delta_time) - 1,
            'ksi_N': 0,
            'ksi_R': 0,
            'r_U': (1 + r_U_daily_target) ** (1. / self.Delta_time) - 1,
            'r_P': 0,
            'r_AP': 0,
            'r_N': (1 + r_N_daily_target) ** (1. / self.Delta_time) - 1,
            'r_R': (1 + r_R_daily_target) ** (1. / self.Delta_time) - 1,
            'd_start_exp': 0.,
            'experiment': "baseline_vaccine_tag"
        }

        self.optimization_general_imperf_testing_with_lockdown = {
            'tau_paramA': (1 + tau_param_A_daily_target) ** (1. / self.Delta_time) - 1,
            'test_sens': 0.75,
            'test_spec': 0.75,
            'ksi_U': 0,
            'ksi_P': (1 + ksi_P_daily_target) ** (1. / self.Delta_time) - 1,
            'ksi_N': 0,
            'ksi_R': 0,
            'r_U': (1 + r_U_daily_target) ** (1. / self.Delta_time) - 1,
            'r_P': 0,
            'r_AP': 0,
            'r_N': (1 + r_N_daily_target) ** (1. / self.Delta_time) - 1,
            'r_R': (1 + r_R_daily_target) ** (1. / self.Delta_time) - 1,
            'd_start_exp': 0.,
            'experiment': "baseline_vaccine_tag"
        }



        self.common_quarantine = {
            'tau_paramA'            : (1+tau_param_A_daily_target)**(1./self.Delta_time)-1,
            'test_sens'     : 1.0,
            'test_spec'     : 1.0,
            'ksi_U'           : (1+ksi_U_daily_target)**(1./self.Delta_time)-1,
            'ksi_P'           : (1+ksi_P_daily_target)**(1./self.Delta_time)-1,
            'ksi_N'           : (1+ksi_N_daily_target)**(1./self.Delta_time)-1,
            'ksi_R'           : (1+ksi_R_daily_target)**(1./self.Delta_time)-1,
            'r_U'           : (1+r_U_daily_target)**(1./self.Delta_time)-1,
            'r_P'           : (1+r_P_daily_target)**(1./self.Delta_time)-1,
            'r_AP'          : (1+r_AP_daily_target)**(1./self.Delta_time)-1,
            'r_N'           : (1+r_N_daily_target)**(1./self.Delta_time)-1,
            'r_R'           : (1+r_R_daily_target)**(1./self.Delta_time)-1,
            'experiment'    : "baseline_vaccine_tag"
        }

    def solve_case(self, model, policy):

        lockdown_policy = policy.lockdown_policy
        testing_policy = policy.testing_policy
        # debugging prints:
        #print("Solving model: ", model)
        #print("Lockdownpolicy: ", lockdown_policy)
        #print("Other params: ", self.__dict__)

        M0_vec = np.zeros(12)
        M0_vec[3] = self.InitialInfect / self.pop  # initial infected, asymptomatic, not quarantined, and unknown cases
        M0_vec[8] = 2. / self.pop  # initial infected, symptomatic, quarantined (and known) cases
        M0_vec[0] = 1 - np.sum(M0_vec)

        Q_inds = [1, 4, 5, 6, 8, 10]
        NQ_inds = [0, 2, 3, 7, 9]
        IANQ_inds = [3] # doesn't include false negatives!
        IAQ_inds = [4, 5]
        ISQ_inds = [8]
        NANQ_inds = [0, 2]
        RANQ_inds = [9]
        NAQ_inds = [1]
        RAQ_inds = [10]

        FP_inds = [6]
        FPQ_inds = [6]
        FN_inds = [7]
        FNNQ_inds = [7]

        test_inds = [0, 1, 3, 4]
        TT_test_inds = [1, 4]

        # Members in each state on different time steps?
        M_t = np.zeros((12, self.T))
        alpha_T = np.zeros((self.T)) # alpha values will be saved in this one
        ksi_TT_T = np.zeros((self.T)) # test and trace Q rate will be saved here
        M_t[:, 0] = M0_vec
        lockdown_effs = np.zeros((self.T))
        tests = np.zeros((self.T))

        def policy_timer(time, policy, default="NA"):
            # policy: dictionary with policy start times as keys
            # time: moment of time at hand
            # returns correct policy parameter value for time

            #debug:
            #print("policy in timer: ", policy)
            if policy == "NA":
                return default
            else:
                filt_keys = []
                for k in list(policy.keys()):
                    if k * self.Delta_time <= time:

                        filt_keys.append(k)

                try:
                    t_key = np.max(filt_keys) # finds the largest key of those <= to time
                    param = policy[t_key]
                except:
                    #print("returning default for policy at time = ", time)
                    return default

            #if time in filt_keys:
                #print("returning param ", param, " for time = ", time)

                return param


        for t in range(1, self.T):

            # calculate policy value for lockdown strength:
            lockdown_eff = policy_timer(t, lockdown_policy, 1.0)

            lockdown_effs[t] = lockdown_eff
            # print("at time ", t, " ksi_U_t = ", ksi_U_t)

            Mt = M_t[:, t - 1]      # compartment 'masses' from last computed time step

            Mt_Q = np.sum(Mt[Q_inds])  # mass of people in quarantine t-1
            Mt_NQ = np.sum(Mt[NQ_inds])  # mass of people out of quarantine t-1

            Mt_IANQ = np.sum(Mt[IANQ_inds])
            Mt_IAQ = np.sum(Mt[IAQ_inds])

            Mt_ISQ = np.sum(Mt[ISQ_inds])

            Mt_NANQ = np.sum(Mt[NANQ_inds])
            Mt_RANQ = np.sum(Mt[RANQ_inds])

            Mt_NAQ = np.sum(Mt[NAQ_inds])
            Mt_RAQ = np.sum(Mt[RAQ_inds])

            Mt_FP = np.sum(Mt[FP_inds])
            Mt_FPQ = np.sum(Mt[FPQ_inds])
            Mt_FN = np.sum(Mt[FN_inds])
            Mt_FNNQ = np.sum(Mt[FNNQ_inds])

            # Masses of groups from which tested persons are selected:
            Mt_test = np.sum(Mt[test_inds])
            Mt_TT_test = np.sum(Mt[TT_test_inds])

            Mt_Total = lockdown_eff*self.lambda_param * Mt_NQ + self.lambda_paramQ * Mt_Q

            Mt_I = lockdown_eff*self.lambda_param * (Mt_IANQ + Mt_FNNQ) + self.lambda_paramQ * (
                        Mt_IAQ + Mt_ISQ)  # added false negatives
            Mt_N = lockdown_eff*self.lambda_param * (Mt_NANQ + Mt_RANQ) + self.lambda_paramQ * (Mt_NAQ + Mt_RAQ + Mt_FPQ)

            # conditional on meeting a person, the probability that they are infected (I) or not (N)
            pit_I = Mt_I / Mt_Total
            pit_N = Mt_N / Mt_Total

            # conditional on meeting an infected person, probability that they are asymptomatic
            pit_IA = (lockdown_eff*self.lambda_param * Mt_IANQ + self.lambda_paramQ * Mt_IAQ + lockdown_eff*self.lambda_param * Mt_FNNQ) / Mt_I  # added false negatives
            pit_IS = (self.lambda_paramQ * Mt_ISQ) / Mt_I

            alphat = pit_I * (pit_IS * self.rhoS + pit_IA * self.rhoA)
            alpha_T[t] = alphat # saves the alpha for this time step


            # A_daily just selects every 14th entry starting at the 14th entry (end of day each day)

            if t <= model['d_start_exp']:
                ksi_U_t = 0
                ksi_P_t = 0
                ksi_N_t = 0
                ksi_R_t = 0

                r_U_t = 0
                r_P_t = 0
                r_AP_t = 0
                r_N_t = 0
                r_R_t = 0

                tau_t = 0
                #tau_re_t = tau_t * self.delta / 2  # TODO: improve: approximation for retesting rate of positives if not symptomatic

                test_sens = 1
                test_spec = 1
                tau_TT = 0

            elif t >= self.d_vaccine:
                ksi_U_t = model['ksi_U']
                ksi_P_t = model['ksi_P']
                ksi_N_t = model['ksi_N']
                ksi_R_t = 0.

                r_U_t = model['r_U']
                r_P_t = model['r_P']
                r_AP_t = model['r_AP']
                r_N_t = model['r_N']
                r_R_t = model['r_R']

                test_sens = model['test_sens']
                test_spec = model['test_spec']
                tau_TT = model['tau_TT']
                #tau_re_t = tau_t * self.delta / 2

            else:
                ksi_U_t = model['ksi_U']
                ksi_P_t = model['ksi_P']
                ksi_N_t = model['ksi_N']
                ksi_R_t = model['ksi_R']

                r_U_t = model['r_U']
                r_P_t = model['r_P']
                r_AP_t = model['r_AP']
                r_N_t = model['r_N']
                r_R_t = model['r_R']

                tau_t = policy_timer(t, testing_policy, model['tau_paramA'])
                tau_TT = model['tau_TT']
                #tau_re_t = tau_t * self.delta / 2
                test_sens = model['test_sens']
                test_spec = model['test_spec']

            # Test and trace rates:
            Mt_tm1 = M_t[:, t - 2]  # compartment 'masses' from previous to last computed time step
            Mtm1_IAQ = np.sum(Mt_tm1[IAQ_inds])
            Mtm1_IANQ = np.sum(Mt_tm1[IANQ_inds])
            Mtm1_NANQ = np.sum(Mt_tm1[NANQ_inds])
            Mtm1_NAQ = np.sum(Mt_tm1[NAQ_inds])

            # TODO: check and confirm lockdown effect in these!

            pit_IST = self.delta*(self.lambda_paramQ * Mtm1_IAQ + lockdown_eff* self.lambda_param * Mtm1_IANQ) / Mt_I
            pit_IAT = (lockdown_eff * self.lambda_param * Mtm1_IANQ + self.lambda_paramQ * Mtm1_IAQ) * tau_t * test_sens / Mt_I
            pit_FP = (lockdown_eff * self.lambda_param * Mtm1_NANQ + self.lambda_paramQ * Mtm1_NAQ) * tau_t * (1 - test_spec) / Mt_N

            eta_I =
            eta_N =

            ksi_TT_I = self.eta * lockdown_eff * self.lambda_param * (pit_I * (pit_IST + pit_IAT) + pit_N * pit_FP) # quarantine probability due to test and trace
            ksi_TT_N =

            ksi_TT_T[t] = ksi_TT

            #print("pit_I:", pit_I)
            #print("pit_IST:", pit_IST)
            #print("pit_IAT:", pit_IAT)
            #print("pit_N:", pit_N)
            #print("pit_FP:", pit_FP)
            #ksi_trace_Q = self.eta * self.lambda_paramQ * (pit_I * (pit_IST + pit_IAT) + pit_N * pit_FP)

            # Create transition matrix and fill it with correct values
            transition_matrix_t = np.zeros((12, 12))

            # from not known NA, NQ - Not infected Asymptomatic, Not Quarantined
            transition_matrix_t[0, 1] = ksi_TT  # To NA, Quarantined, based on test & trace
            transition_matrix_t[0, 2] = tau_t * test_spec  # To known not-infected asymptomatic, NQ
            transition_matrix_t[0, 3] = lockdown_eff * self.lambda_param * alphat  # To unknown infected asympt., not NQ
            transition_matrix_t[0, 6] = tau_t * (1.0 - test_spec)  # to false positive, NQ


            # from not known NA, Q - Not infected Asymptomatic, Quarantined
            transition_matrix_t[1, 0] = r_U_t
            transition_matrix_t[1, 2] = (tau_t + tau_TT) * test_spec  # To known not-infected asymptomatic, Quarantined
            transition_matrix_t[1, 4] = self.lambda_paramQ * alphat
            transition_matrix_t[1, 6] = (tau_t + tau_TT) * (1.0 - test_spec)  # To false positive, Quarantined


            # from known NA, NQ - Not infected Asymptomatic, Not Quarantined
            transition_matrix_t[2, 0] = self.sigma
            transition_matrix_t[2, 3] = lockdown_eff * self.lambda_param * alphat  # To unknown infected asymptomatic, not NQ


            # from not known IA, NQ - Infected Asymptomatic, Not Quarantined
            transition_matrix_t[3, 4] = ksi_TT
            transition_matrix_t[3, 5] = tau_t * test_sens  # To known infected asymptomatic, Q
            transition_matrix_t[3, 7] = tau_t * (1.0 - test_sens)  # To false negative, NQ
            transition_matrix_t[3, 8] = self.delta
            transition_matrix_t[3, 9] = self.omegaR

            # from not known IA, Q - Infected Asymptomatic, Quarantined
            transition_matrix_t[4, 3] = r_U_t
            transition_matrix_t[4, 5] = tau_TT * test_sens  # To known infected asymptomatic, Quarantined
            transition_matrix_t[4, 7] = tau_TT * (1.0 - test_sens)  # to false negative, Quarantined
            transition_matrix_t[4, 8] = self.delta
            transition_matrix_t[4, 10] = self.omegaR

            # from known IA, Q - Infected Asymptomatic, Quarantined
            transition_matrix_t[5, 8] = self.delta
            transition_matrix_t[5, 10] = self.omegaR

            # from false Positive, Quarantined (index 9)
            # i.e. not infected asymptomatic but treated like infected
            transition_matrix_t[6, 0] = self.omegaR
            transition_matrix_t[6, 5] = self.lambda_paramQ * alphat  # to known infected, quarantined (actually gets infected) 'infection while in Q rate'

            # from False Negative, Not Quarantined (index 10)
            # i.e. infected (asymptomatic) but treated like not infected

            transition_matrix_t[7, 8] = self.delta  # to infected symptomatic not quarantined  - assume infection diagnosed correctly then - "symptom dev rate"
            transition_matrix_t[7, 9] = self.omegaR


            # from (known) Infected Symptomatic, Quarantined
            transition_matrix_t[8, 10] = self.omegaR
            transition_matrix_t[8, 11] = self.omegaD

            # from Recovered Asymptomatic, Not Quarantined
            transition_matrix_t[9, 0] = self.gamma  # immunity loss
            transition_matrix_t[9, 10] = ksi_TT

            # from Recovered Asymptomatic, Quarantined
            transition_matrix_t[10, 9] = r_R_t

            # probabilities for staying in same compartment
            transition_matrix_t += np.diag(1 - np.sum(transition_matrix_t, axis=1))

            # This tests that there are no clearly faulty values in the matrix

            assert np.min(transition_matrix_t) >= 0
            assert np.max(transition_matrix_t) <= 1

            # M_t at t calculated from previous time step Mt and transitions thru matrix multiplication
            # .T is transpose
            # @ is matrix multiplication
            M_t[:, t] = transition_matrix_t.T @ Mt
            tests[t] = (Mt_test*tau_t + Mt_TT_test*tau_TT)*self.pop

        # Total productivity = productivity of non quarantined + productivity of quarantined non-symptomatic
        Y_t = lockdown_effs * np.sum(M_t[[0, 2, 3, 7, 9]], axis=0) + \
              self.A_rel * np.sum(M_t[[1, 4, 5, 6, 10]], axis=0)

        Reported_T_start = self.pop * ((test_sens * tau_t + self.delta) * M_t[3] + (test_sens * tau_TT + self.delta) * M_t[4] + test_spec * tau_t * M_t[0] + test_spec * (tau_t + tau_TT) * M_t[1] ) # reported calculated from tested cases
        Reported_T_start[0] = 0
        Reported_T = np.cumsum(Reported_T_start)

        Reported_D = Reported_T[13::14]  # Note: 13::14 refers to time indices, i.e. 'end of day for all days'
        Notinfected_D = np.sum(M_t[[0, 1, 2]], axis=0)[13::14]
        Unreported_D = np.sum(M_t[[3, 4]], axis=0)[13::14]
        Infected_D = (np.sum(M_t[[3, 4, 5]], axis=0) + np.sum(M_t[[7, 8]], axis=0))[13::14]
        Infected_in_Q = np.sum(M_t[[4, 5, 8]], axis=0)[
                        13::14]  # includes all infected in quarantine including false negs
        Infected_not_Q = np.sum(M_t[[3, 7]], axis=0)[13::14]  # includes false negatives
        Symptomatic_D = np.sum(M_t[ISQ_inds], axis=0)[13::14]
        False_pos = np.sum(M_t[[6]], axis=0)[13::14]
        False_neg = np.sum(M_t[[7]], axis=0)[13::14]
        Recovered_D = np.sum(M_t[[9,10]], axis=0)[13::14]
        Dead_D = M_t[11][13::14]    # Dead at end of each day
        Infected_T = np.sum(M_t[[3, 4, 5]], axis=0) + np.sum(M_t[[7, 8]], axis=0)
        Y_D = Y_t[13::14]
        Y_total = np.sum(Y_t)
        total_cost = sum(tests)*self.test_cost

        Unk_NA_nQ_D = M_t[0][13::14]
        Unk_NA_Q_D = M_t[1][13::14]
        K_NA_nQ_D = M_t[2][13::14]
        Unk_IA_nQ_D = M_t[3][13::14]
        Unk_IA_Q_D = M_t[4][13::14]
        K_IA_Q_D = M_t[5][13::14]
        tests_D = tests[13::14]
        ksi_TT_D = ksi_TT_T[13::14]
        alpha_D = alpha_T[13::14]

        return Reported_D, Notinfected_D, Unreported_D, Infected_D, \
               False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total, total_cost, tests_D, Unk_NA_nQ_D, Unk_NA_Q_D, K_NA_nQ_D, Unk_IA_nQ_D, Unk_IA_Q_D, K_IA_Q_D, alpha_D, ksi_TT_D, Symptomatic_D

    def solve_model(self, lockdown_policy={10000: 0}, testing_policy = {10000: 0}):
        Reported_D_base, Notinfected_D_base, Unreported_D_base, Infected_D_base, \
                False_pos_base, False_neg_base, Recovered_D_base, Dead_D_base, Infected_T_base, Infected_not_Q_base, Infected_in_Q_base, Y_D_base, M_t_base, Y_total_base, total_cost_base = \
                self.solve_case(self.baseline, lockdown_policy, testing_policy)
        Tstar = np.argwhere(Reported_D_base>100)[0][0]
        YearsPlot = 3
        Tplot = np.arange(Tstar, min(Tstar + YearsPlot * 365, self.T/self.Delta_time) + .5, 1)
        Xplot = np.arange(0, len(Tplot))
        self.Tstar = Tstar

        self.common_quarantine['d_start_exp'] = (Tstar+1) * self.Delta_time + \
                self.policy_offset * self.Delta_time

        Reported_D_com, Notinfected_D_com, Unreported_D_com, Infected_D_com, \
                False_pos_com, False_neg_com, Recovered_D_com, Dead_D_com, Infected_T_com, Infected_not_Q_com, Infected_in_Q_com, Y_D_com, M_t_com, Y_total_com, total_cost_com = \
                self.solve_case(self.common_quarantine, lockdown_policy, testing_policy)

        return Reported_D_com, Infected_D_com, Dead_D_com, Y_D_com, False_pos_com, False_neg_com, Infected_not_Q_com, Infected_in_Q_com, Y_total_com, total_cost_com

    def run_experiment(self, tau_param, Delta, test_sens, test_spec, lockdown_policy={10000: 0.0}, testing_policy={10000: 0.0}):

        tau_param_A_daily_target = tau_param

        r_U_daily_target	= 0
        r_N_daily_target	= 0
        r_P_daily_target	= 0
        r_R_daily_target	= self.r_high
        r_AP_daily_target   = self.r_AP # self.r

        ksi_U_daily_target   = self.ksi_base
        ksi_P_daily_target   = self.ksi_base_high
        ksi_N_daily_target   = self.ksi_base*Delta
        ksi_R_daily_target   = 0
        test_sens_exp      = test_sens
        test_spec_exp      = test_spec

        self.test_and_quarantine = {
            'tau_paramA'            : (1+tau_param_A_daily_target)**(1./self.Delta_time)-1,
            'test_sens'     : test_sens_exp,
            'test_spec'     : test_spec_exp,
            'ksi_U'           : (1+ksi_U_daily_target)**(1./self.Delta_time)-1,
            'ksi_P'           : (1+ksi_P_daily_target)**(1./self.Delta_time)-1,
            'ksi_N'           : (1+ksi_N_daily_target)**(1./self.Delta_time)-1,
            'ksi_R'           : (1+ksi_R_daily_target)**(1./self.Delta_time)-1,
            'r_U'           : (1+r_U_daily_target)**(1./self.Delta_time)-1,
            'r_P'           : (1+r_P_daily_target)**(1./self.Delta_time)-1,
            'r_AP'          : (1 + r_AP_daily_target) ** (1. / self.Delta_time) - 1,
            'r_N'           : (1+r_N_daily_target)**(1./self.Delta_time)-1,
            'r_R'           : (1+r_R_daily_target)**(1./self.Delta_time)-1,
            'experiment'    : "baseline_vaccine_tag"
        }

        self.test_and_quarantine['d_start_exp'] = (self.Tstar+1) * self.Delta_time + \
                self.policy_offset * self.Delta_time

        # NOTE: here base case used instead of original test_and_quarantine
        Reported_D_test, Notinfected_D_test, Unreported_D_test, Infected_D_test, \
                False_pos_test, False_neg_test, Recovered_D_test, Dead_D_test, Infected_T_test,  Infected_not_Q_test, Infected_in_Q_test, Y_D_test, M_t_test, Y_total_test, total_cost_test = \
                self.solve_case(self.optimization_no_test, lockdown_policy, testing_policy)

        return Reported_D_test, Infected_D_test, Dead_D_test, Y_D_test, False_pos_test, False_neg_test, Infected_not_Q_test, Infected_in_Q_test, Y_total_test, total_cost_test



def generate_plots(Delta, tau_param, test_sens, test_spec, ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
             omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect, slide_var, lockdown_policy, testing_policy):

    if slide_var == 6:  # slide over lockdown policies
        active_policies = testing_policy

    elif slide_var == 5:
        active_policies = lockdown_policy


    colors = ['red', 'blue']
    styles = ['dot', 'dash']

    rmin = 0
    rmax = 0
    imin = 0
    imax = 0
    dmin = 0
    dmax = 0
    ymin = .5
    ymax = 0
    fpmin = 0
    fpmax = 0
    fnmin = 0
    fnmax = 0
    inqmin = 0
    inqmax = 0
    iqmin = 0
    iqmax = 0

    fig = make_subplots(5, 2, print_grid = False, \
                        subplot_titles=("A. Reported cases", "B. Current symptomatic cases", "C. Deaths - Cumulative", "D. Current output", "E. False positives", "F: False negatives", "G: Infected, not quarantined", "H: Infected, in quarantine", "I: Lockdown policy", "J: optimization outcomes"),
                        vertical_spacing = .2, specs=[[{},{}], [{},{}], [{},{}], [{},{}], [{}, {"secondary_y": True}]])

    print("Creating a corona model with testing rate = ", tau_param, " sensitivity = ", test_sens, " and specificity = ", test_spec)
    model = optimizable_corona_model(ksi_base, A_rel, r_AP, d_vaccine, rel_rho, delta_param, \
                 omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect)

    Reported_D_com, Infected_D_com, Dead_D_com, Y_D_com, False_pos_com, False_neg_com, Infected_not_Q_com, Infected_in_Q_com, Y_total_com, total_cost_com = model.solve_model()

    rmin = min(rmin, np.min(Reported_D_com) * 1.2)
    rmax = max(rmax, np.max(Reported_D_com) * 1.2)
    imin = min(imin, np.min(Infected_D_com) * 1.2)
    imax = max(imax, np.max(Infected_D_com) * 1.2)
    dmin = min(dmin, np.min(Dead_D_com) * 1.2)
    dmax = max(dmax, np.max(Dead_D_com) * 1.2)
    ymin = min(ymin, np.min(Y_D_com) * 1.2)
    ymax = max(ymax, np.max(Y_D_com) * 1.2)
    fpmin = min(fpmin, np.min(False_pos_com) * 1.2)
    fpmax = max(fpmax, np.max(False_pos_com) * 1.2)
    fnmin = min(fnmin, np.min(False_neg_com) * 1.2)
    fnmax = max(fnmax, np.max(False_neg_com) * 1.2)
    inqmin = min(inqmin, np.min(Infected_not_Q_com) * 1.2)
    inqmax = max(inqmax, np.max(Infected_not_Q_com) * 1.2)
    iqmin = min(iqmin, np.min(Infected_in_Q_com) * 1.2)
    iqmax = max(iqmax, np.max(Infected_in_Q_com) * 1.2)
    pmin = 0
    pmax = np.max(list(active_policies[0].values()))


    outmin = 0.0
    outmax = Y_total_com

    xticks = list(active_policies[0].keys())

    fig.add_scatter(y = Reported_D_com, row = 1, col = 1, visible = True, showlegend = True,
                    name = 'Base case', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Infected_D_com, row = 1, col = 2, visible = True, showlegend = False,
                    name = 'Base case', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Dead_D_com, row = 2, col = 1, visible = True, showlegend = False,
                    name = 'Base case', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Y_D_com, row = 2, col = 2, visible = True, showlegend = False,
                    name = 'Base case', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y=False_pos_com, row=3, col=1, visible=True, showlegend=False,
                    name='Base case', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=False_neg_com, row=3, col=2, visible=True, showlegend=False,
                    name='Base case', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=Infected_not_Q_com, row=4, col=1, visible=True, showlegend=False,
                    name='Base case', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=Infected_in_Q_com, row=4, col=2, visible=True, showlegend=False,
                    name='Base case', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_trace(go.Bar(x = ['Output, baseline'], y = [Y_total_com], width=[0.5]), row=5, col=2,
                  secondary_y=True),
    fig.add_trace(go.Bar(x=['Deaths, baseline'], y=[Dead_D_com[-1]], width=[0.5]), row=5, col=2, secondary_y=False),
    #fig.add_trace(go.Bar(y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), row=5, col=1)
    fig.add_trace(go.Bar(x=xticks, y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  row=5, col=1)




    if slide_var == 1: #Slide over tau_param
        prd = product(tau_param, [Delta], [test_sens], [test_spec], [lockdown_policy], [testing_policy])
        slider_vars = tau_param
        slider_varname = "tau_param"

    if slide_var == 2: #Slide over Delta
        prd = product([tau_param], Delta, [test_sens], [test_spec], [lockdown_policy], [testing_policy])
        slider_vars = Delta
        slider_varname = "Delta"

    if slide_var == 3:  # Slide over test_sens
        prd = product([tau_param], [Delta], test_sens, [test_spec], [lockdown_policy], [testing_policy])
        slider_vars = test_sens
        slider_varname = "test sensitivity"

    if slide_var == 4:  # Slide over test_spec
        prd = product([tau_param], [Delta], [test_sens], test_spec, [lockdown_policy], [testing_policy])
        slider_vars = test_spec
        slider_varname = "test specificity"

    if slide_var == 5: # slide over lockdown policies
        prd = product([tau_param], [Delta], [test_sens], [test_spec], lockdown_policy, testing_policy)
        slider_vars = range(0,len(lockdown_policy))
        slider_varname = "lockdown policy"
        active_policies = lockdown_policy

    if slide_var == 6: # slide over lockdown policies
        prd = product([tau_param], [Delta], [test_sens], [test_spec], lockdown_policy, testing_policy)
        slider_vars = range(0,len(testing_policy))
        slider_varname = "testing policy"
        active_policies = testing_policy


    pool = Pool(os.cpu_count())
    print("starting experiment")
    results = pool.starmap(model.run_experiment, prd)


    for j in range(len(slider_vars)):
        print("Output from results: ", results[j][8])

        rmin = min(rmin, np.min(results[j][0]) * 1.2)
        rmax = max(rmax, np.max(results[j][0]) * 1.2)
        imin = min(imin, np.min(results[j][1]) * 1.2)
        imax = max(imax, np.max(results[j][1]) * 1.2)
        dmin = min(dmin, np.min(results[j][2]) * 1.2)
        dmax = max(dmax, np.max(results[j][2]) * 1.2)
        ymin = min(ymin, np.min(results[j][3]) * 1.2)
        ymax = max(ymax, np.max(results[j][3]) * 1.2)
        fpmin = min(fpmin, np.min(results[j][4]) * 1.2)
        fpmax = max(fpmax, np.max(results[j][4]) * 1.2)
        fnmin = min(fnmin, np.min(results[j][5]) * 1.2)
        fnmax = max(fnmax, np.max(results[j][5]) * 1.2)
        inqmin = min(inqmin, np.min(results[j][6]) * 1.2)
        inqmax = max(inqmax, np.max(results[j][6]) * 1.2)
        iqmin = min(iqmin, np.min(results[j][7]) * 1.2)
        iqmax = max(iqmax, np.max(results[j][7]) * 1.2)
        pmax = max(pmax, np.max(list(active_policies[j].values())) * 1.2)
        outmax = max(outmax, results[j][8])



        fig.add_scatter(y = results[j][0], row = 1, col = 1, visible = j == 0, showlegend = True,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash = styles[1]))
        fig.add_scatter(y = results[j][1], row = 1, col = 2, visible = j == 0, showlegend = False,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash = styles[1]))
        fig.add_scatter(y = results[j][2], row = 2, col = 1, visible = j == 0, showlegend = False,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash = styles[1]))
        fig.add_scatter(y = results[j][3], row = 2, col = 2, visible = j == 0, showlegend = False,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash = styles[1]))
        fig.add_scatter(y = results[j][4], row=3, col = 1, visible = j == 0, showlegend=False,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash=styles[1]))
        fig.add_scatter(y = results[j][5], row=3, col = 2, visible = j == 0, showlegend=False,
                        name = 'Quarantine & Test', line = dict(color = (colors[1]), width = 3, dash=styles[1]))
        fig.add_scatter(y=results[j][6], row=4, col=1, visible=j == 0, showlegend=False,
                        name='Quarantine & Test', line=dict(color=(colors[1]), width=3, dash=styles[1]))
        fig.add_scatter(y=results[j][7], row=4, col=2, visible=j == 0, showlegend=False,
                        name='Quarantine & Test', line=dict(color=(colors[1]), width=3, dash=styles[1]))
        #fig.add_trace(go.Bar(y=list(lockdown_policy[j].values())), row=5, col=1)

        fig.add_trace(go.Bar(x = ['Output'], y=[results[j][8]], width=[0.5]), row=5, col=2, secondary_y=True),
        fig.add_trace(go.Bar(x=['Deaths'], y=[results[j][2][-1]], width=[0.5]),
                      row=5, col=2, secondary_y=False),
        fig.add_trace(go.Bar(x=xticks, y=list(active_policies[j].values())),
                      row=5, col=1)

    steps = []
    n_plots = 11 # number of subplots in fig. NOTE: secondary axes also count as a new subplot
    for i in range(len(slider_vars)):
        step = dict(
            method = 'restyle',
            args = [{'visible': ['legendonly'] * len(fig.data)},
                    {'showlegend': ['False'] * len(fig.data)}],
            label = slider_varname + ' = \n'+'{}'.format(round(slider_vars[i], 3))
        )
        step['args'][1]['showlegend'][0] = True
        step['args'][1]['showlegend'][n_plots + i * n_plots] = True
        for j in range(n_plots):
            step['args'][0]['visible'][int(j)] = True
        for j in range(n_plots):
            try:
                step['args'][0]['visible'][n_plots + j + i * n_plots] = True
            except IndexError:
                print("skipped index: ", n_plots + j + i * n_plots)

        steps.append(step)

    sliders = [dict(
        steps = steps
    )]

    fig.layout.sliders = sliders
    for i in fig['layout']['annotations']:
        i['font'] = dict(color='black', size = 16)
    fig['layout'].update(height=1000, width=1000, showlegend = False)

    fig['layout']['xaxis1'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 365*model.T_years], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis2'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 365*model.T_years], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis3'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 365*model.T_years], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis4'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 365*model.T_years], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis5'].update(title=go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 365*model.T_years], \
                                    gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis6'].update(title=go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 365*model.T_years], \
                                    gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis7'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 365*model.T_years], \
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis8'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 365*model.T_years], \
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)

    fig['layout']['xaxis9'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 365*model.T_years], \
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True )
    fig['layout']['xaxis10'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 4], \
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)



    fig['layout']['yaxis1'].update(title=go.layout.yaxis.Title(
                                text='Logarithm - Base 10', font=dict(color='black')), type='log', range = [rmin, np.log10(rmax)], gridcolor = 'rgb(220,220,220)', \
                                   showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis2'].update(title=go.layout.yaxis.Title(
                                text='Fraction of Initial Population', font=dict(color='black')), range=[imin, imax], gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis3'].update(title=go.layout.yaxis.Title(
                                text='Fraction of Initial Population', font=dict(color='black')), range = [dmin, dmax], gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis4'].update(title=go.layout.yaxis.Title(
                                text='Output', font=dict(color='black')), range = [ymin, 1.05], gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis5'].update(title=go.layout.yaxis.Title(
                                text='Fraction of Initial Population', font=dict(color='black')), range=[fpmin, fpmax], gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis6'].update(title=go.layout.yaxis.Title(
                                text='Fraction of Initial Population', font=dict(color='black')), range=[fnmin, fnmax], gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis7'].update(title=go.layout.yaxis.Title(
        text='Fraction of Initial Population', font=dict(color='black')), range=[inqmin, inqmax],
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis8'].update(title=go.layout.yaxis.Title(
        text='Fraction of Initial Population', font=dict(color='black')), range=[iqmin, iqmax],
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)

    fig['layout']['yaxis9'].update(title=go.layout.yaxis.Title(
        text='Parameter value for xi^U', font=dict(color='black')), range=[pmin, pmax],
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis10'].update(title=go.layout.yaxis.Title(
        text='Deaths', font=dict(color='black')), range=[dmin, dmax],
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['yaxis11'].update(title=go.layout.yaxis.Title(
        text='Output', font=dict(color='black')), range=[outmin, outmax],
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)

    fig['layout']['margin'].update(l=20, r=20, t=20, b=20)

    fig['layout']['plot_bgcolor'] = 'rgba(0,0,0,0)'

    return fig

def generate_plots_2d(Delta, tau_param, test_sens, test_spec, ksi_base, A_rel, d_vaccine, rel_rho, delta_param, \
             omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect):

    model = optimizable_corona_model(ksi_base, A_rel, d_vaccine, rel_rho, delta_param, \
                 omegaR_param, pii_D, R_0, rel_lambda_param, initial_infect)

    Reported_D, Infected_D, Dead_D, Y_D, Reported_D_com, Infected_D_com, \
        Dead_D_com, Y_D_com, Y_total_com = model.solve_model()

    prd = product(tau_param, Delta)

    pool = Pool(os.cpu_count())
    results = pool.starmap(model.run_experiment, prd)

    return Reported_D, Infected_D, Dead_D, Y_D, Reported_D_com, Infected_D_com, \
        Dead_D_com, Y_D_com, results, prd


