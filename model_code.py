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

class corona_model(object):
    # ξ_base: baseline quarantine rate
    # A_rel: relative productivity of quarantined
    # d_vaccine: date of vaccine
    # rel_ρ: relative infectiousness?
    # δ_param: 1/ mean days to show symptoms
    # ωR_param ~mean days to recovery (check details)
    # π_D: mortality rate
    # R_0: basic reproduction number
    # rel_λ: effectiveness of quarantine?
    # initial_infect

    def __init__(self, ξ_base, A_rel, r_AP, d_vaccine, rel_ρ, δ_param, \
                 ωR_param, π_D, R_0, rel_λ,initial_infect):
        self.pop        = 340_000_000
        self.T_years    = 3
        self.Δ_time     = 14
        self.T          = self.T_years * 365 * self.Δ_time
        self.λ          = 1
        self.γ          = 0         # infection rate multiplier for recovered patients' reinfection
        self.ξ_base_high= .999
        self.r_high     = .999
        self.r          = .98
        self.r_AP       = r_AP
        self.δ          = 1/(self.Δ_time*δ_param)
        self.ωR         = 1/(self.Δ_time*ωR_param)

        self.ωD         = self.ωR*π_D/(1-π_D)
        self.λQ         = rel_λ*self.λ
        self.ρS         = R_0/((self.λ/self.δ)*(rel_ρ + self.δ/(self.ωR+self.ωD)))
        self.ρA         = rel_ρ*self.ρS

        self.InitialInfect = initial_infect
        self.d_vaccine     = d_vaccine
        self.A_rel         = A_rel
        self.ξ_base        = ξ_base

        #self.test_sens      = test_sens
        #self.test_spec      = test_spec

        self.baseline = {
            'τA'            : 0.,
            'test_sens'     : 1.0,
            'test_spec'     : 1.0,
            'ξ_U'           : 0., # baseline quarantine rate
            'ξ_P'           : 0.,
            'ξ_N'           : 0.,
            'ξ_R'           : 0.,
            'r_U'           : self.r,
            'r_P'           : 0.,
            'r_AP'          : self.r_AP,
            'r_N'           : self.r,
            'r_R'           : self.r_high,
            'd_start_exp'   : 0.,
            'experiment'    : "baseline_vaccine_tag"
        }

        τ_A_daily_target    = 0
        ξ_U_daily_target	= ξ_base
        ξ_P_daily_target	= self.ξ_base_high
        ξ_N_daily_target	= ξ_base
        ξ_R_daily_target	= 0

        r_U_daily_target	= 0
        r_N_daily_target	= 0
        r_P_daily_target	= 0
        r_AP_daily_target   = self.r_AP # self.r
        r_R_daily_target	= self.r_high

        self.policy_offset = 14

        self.common_quarantine = {
            'τA'            : (1+τ_A_daily_target)**(1./self.Δ_time)-1,
            'test_sens'     : 1.0,
            'test_spec'     : 1.0,
            'ξ_U'           : (1+ξ_U_daily_target)**(1./self.Δ_time)-1,
            'ξ_P'           : (1+ξ_P_daily_target)**(1./self.Δ_time)-1,
            'ξ_N'           : (1+ξ_N_daily_target)**(1./self.Δ_time)-1,
            'ξ_R'           : (1+ξ_R_daily_target)**(1./self.Δ_time)-1,
            'r_U'           : (1+r_U_daily_target)**(1./self.Δ_time)-1,
            'r_P'           : (1+r_P_daily_target)**(1./self.Δ_time)-1,
            'r_AP'          : (1+r_AP_daily_target)**(1./self.Δ_time)-1,
            'r_N'           : (1+r_N_daily_target)**(1./self.Δ_time)-1,
            'r_R'           : (1+r_R_daily_target)**(1./self.Δ_time)-1,
            'experiment'    : "baseline_vaccine_tag"
        }

    def solve_case(self, model):
        M0_vec = np.zeros(17)
        M0_vec[4] = self.InitialInfect / self.pop # initial infected, asymptomatic, not quarantined, and unknown cases
        M0_vec[8] = 1. / self.pop # initial infected, symptomatic, not quarantined (and known) cases
        M0_vec[0] = 1 - np.sum(M0_vec)

        Q_inds      = [1,3,5,7,9,11,13,15]
        NQ_inds     = [0,2,4,6,8,10,12,14]
        IANQ_inds   = [4,6]
        IAQ_inds    = [5,7]
        ISNQ_inds   = [12]
        ISQ_inds    = [13]
        NANQ_inds   = [0,2]
        RANQ_inds   = [14]
        NAQ_inds    = [1,3]
        RAQ_inds    = [15]

        FP_inds = [8,9]
        FPNQ_inds = [8]
        FPQ_inds = [9]
        FN_inds = [10,11]
        FNNQ_inds = [10]
        FNQ_inds = [11]

        # Members in each state on different time steps?
        M_t = np.zeros((17, self.T))
        M_t[:,0] = M0_vec

        for t in range(1,self.T):
            Mt = M_t[:,t-1]

            Mt_Q        = np.sum(Mt[Q_inds])    # mass of people in quarantine t-1
            Mt_NQ       = np.sum(Mt[NQ_inds])   # mass of people out of quarantine t-1

            Mt_IANQ     = np.sum(Mt[IANQ_inds])
            Mt_IAQ      = np.sum(Mt[IAQ_inds])

            Mt_ISNQ     = np.sum(Mt[ISNQ_inds])
            Mt_ISQ      = np.sum(Mt[ISQ_inds])

            Mt_NANQ     = np.sum(Mt[NANQ_inds])
            Mt_RANQ     = np.sum(Mt[RANQ_inds])

            Mt_NAQ      = np.sum(Mt[NAQ_inds])
            Mt_RAQ      = np.sum(Mt[RAQ_inds])

            Mt_FP       = np.sum(Mt[FP_inds])
            Mt_FPNQ = np.sum(Mt[FPNQ_inds])
            Mt_FPQ = np.sum(Mt[FPQ_inds])
            Mt_FN = np.sum(Mt[FN_inds])
            Mt_FNNQ = np.sum(Mt[FNNQ_inds])
            Mt_FNQ = np.sum(Mt[FNQ_inds])

            Mt_Total    = self.λ*Mt_NQ + self.λQ*Mt_Q
            Mt_I        = self.λ*(Mt_IANQ + Mt_ISNQ + Mt_FNNQ) + self.λQ*(Mt_IAQ + Mt_ISQ + Mt_FNQ) # added false negatives
            Mt_N        = self.λ*(Mt_NANQ + Mt_RANQ + Mt_FPNQ) + self.λQ*(Mt_NAQ + Mt_RAQ + Mt_FPQ)

            pit_I       = Mt_I/Mt_Total
            pit_IA      = (self.λ*Mt_IANQ + self.λQ*Mt_IAQ + self.λ*Mt_FNNQ + self.λQ*Mt_FNQ)/Mt_I # added false negatives
            pit_IS      = (self.λ*Mt_ISNQ + self.λQ*Mt_ISQ)/Mt_I

            alphat      = pit_I*(pit_IS*self.ρS + pit_IA*self.ρA)

            # A_daily just selects every 14th entry starting at the 14th entry (end of day each day)

            if t <= model['d_start_exp']:
                ξ_U_t = 0
                ξ_P_t = 0
                ξ_N_t = 0
                ξ_R_t = 0

                r_U_t = 0
                r_P_t = 0
                r_AP_t = 0
                r_N_t = 0
                r_R_t = 0
                
                tau_t = 0
                test_sens = 1
                test_spec = 1

            elif t >= self.d_vaccine:
                ξ_U_t = model['ξ_U']
                ξ_P_t = model['ξ_P']
                ξ_N_t = model['ξ_N']
                ξ_R_t = 0.

                r_U_t = model['r_U']
                r_P_t = model['r_P']
                r_AP_t = model['r_AP']
                r_N_t = model['r_N']
                r_R_t = model['r_R']

                test_sens = model['test_sens']
                test_spec = model['test_spec']

            else:
                ξ_U_t = model['ξ_U']
                ξ_P_t = model['ξ_P']
                ξ_N_t = model['ξ_N']
                ξ_R_t = model['ξ_R']

                r_U_t = model['r_U']
                r_P_t = model['r_P']
                r_AP_t = model['r_AP']
                r_N_t = model['r_N']
                r_R_t = model['r_R']

                tau_t = model['τA']
                test_sens = model['test_sens']
                test_spec = model['test_spec']

            # Create transition matrix and fill it with correct values
            transition_matrix_t         = np.zeros((17,17))

            # from not known NA, NQ - Not infected Asymptomatic, Not Quarantined
            transition_matrix_t[0,1]    = ξ_U_t                     # To NA, Quarantined
            transition_matrix_t[0,2]    = tau_t * test_spec         # To known not-infected asymptomatic, NQ
            transition_matrix_t[0,8]    = tau_t * (1.0 - test_spec)   # to false positive, NQ
            transition_matrix_t[0,4]    = self.λ*alphat             # To unknown infected asymptomatic, not NQ


            # from not known NA, NQ - Not infected Asymptomatic, Quarantined
            transition_matrix_t[1,0]    = r_U_t
            transition_matrix_t[1,3]    = tau_t*test_spec           # To known not-infected asymptomatic, Quarantined
            transition_matrix_t[1,9]    = tau_t*(1.0-test_spec)       # To false positive, Quarantined
            transition_matrix_t[1,5]    = self.λQ*alphat

            # from known NA, NQ - Not infected Asymptomatic, Not Quarantined
            transition_matrix_t[2,3]    = ξ_N_t
            transition_matrix_t[2,6]    = self.λ*alphat

            # from known NA, NQ - Not infected Asymptomatic, Quarantined
            transition_matrix_t[3,2]    = r_N_t
            transition_matrix_t[3,7]    = self.λQ*alphat

            # from not known IA, NQ - Infected Asymptomatic, Not Quarantined
            transition_matrix_t[4,5]    = ξ_U_t
            transition_matrix_t[4,6]    = tau_t*test_sens           # To known infected asymptomatic, NQ
            transition_matrix_t[4,10]   = tau_t*(1.0-test_sens)       # To false negative, NQ
            transition_matrix_t[4,12]    = self.δ

            # from not known IA, Q - Infected Asymptomatic, Quarantined
            transition_matrix_t[5,4]    = r_U_t
            transition_matrix_t[5,7]    = tau_t*test_sens           # To known infected asymptomatic, Quarantined
            transition_matrix_t[5,11]   = tau_t*(1.0-test_sens)       # to false negative, Quarantined
            transition_matrix_t[5,13]    = self.δ

            # from known IA, NQ - Infected Asymptomatic, Not Quarantined
            transition_matrix_t[6,7]    = ξ_P_t
            transition_matrix_t[6,12]    = self.δ

            # from known IA, Q - Infected Asymptomatic, Quarantined
            transition_matrix_t[7,6]    = r_AP_t
            transition_matrix_t[7,13]    = self.δ

            # from false Positive, Not Quarantined (index 8)
            # i.e. not infected asymptomatic but treated like infected

            transition_matrix_t[8, 6] = self.λ * alphat  # to known infected, not quarantined (actually gets infected) - 'infection while not in Q rate'
            transition_matrix_t[8, 9] = ξ_P_t  # to false positive, quarantined - 'quarantine rate for known infected'

            # from false Positive, Quarantined (index 9)
            # i.e. not infected asymptomatic but treated like infected

            transition_matrix_t[9, 7] = self.λQ * alphat  # to known infected, quarantined (actually gets infected) 'infection while in Q rate'
            transition_matrix_t[9, 8] = r_AP_t  # to false positive not quarantined  - 'quarantine release rate'

            # from False Negative, Not Quarantined (index 10)
            # i.e. infected (asymptomatic) but treated like not infected

            transition_matrix_t[10, 12] = self.δ  # to infected symptomatic not quarantined  - assume infection diagnosed correctly then - "symptom dev rate"
            transition_matrix_t[10, 11] = ξ_N_t  # to false negative, quarantined - 'known not infected quarantine rate'
            # transition_matrix_t[10, 14] = self.ωR  # to recovered, not quarantined
            # transition_matrix_t[10, 16] = self.ωD  # death due COVID-19

            # from False Negative, Quarantined (index 11)
            # i.e. infected (asymptomatic) but treated like not infected, but quarantined

            transition_matrix_t[11, 13] = self.δ  # to infected symptomatic quarantined - assume infection diagnosed correctly then?
            transition_matrix_t[11, 10] = r_N_t  # to false negative, not quarantined - 'quarantine release rate'
            # transition_matrix_t[11, 15] = self.ωR  # to recovered, quarantined
            # transition_matrix_t[11, 16] = self.ωD  # death due COVID-19

            # from (known) Infected Symptomatic, Not Quarantined
            transition_matrix_t[12,13]    = ξ_P_t
            transition_matrix_t[12,14]   = self.ωR
            transition_matrix_t[12,16]   = self.ωD

            # from (known) Infected Symptomatic, Quarantined
            transition_matrix_t[13,12]    = r_P_t
            transition_matrix_t[13,15]   = self.ωR
            transition_matrix_t[13,16]   = self.ωD


            # from Recovered Asymptomatic, Not Quarantined
            transition_matrix_t[14,4]    = self.γ * self.λ*alphat   # reinfection - so far has been 0
            transition_matrix_t[14,15]   = ξ_R_t

            # from Recovered Asymptomatic, Quarantined
            transition_matrix_t[15,5]    = self.γ * self.λ*alphat   # reinfection - so far has been 0
            transition_matrix_t[15,14]   = r_R_t

            if t >= self.d_vaccine:
                transition_matrix_t[0,14] = .001
                transition_matrix_t[1,14] = .001
                transition_matrix_t[2,14] = .001
                transition_matrix_t[3,14] = .001
                transition_matrix_t[8,14] = .001
                transition_matrix_t[9,14] = .001

            transition_matrix_t += np.diag(1 - np.sum(transition_matrix_t, axis=1))

            # This tests that there are no clearly faulty values in the matrix

            assert np.min(transition_matrix_t) >= 0
            assert np.max(transition_matrix_t) <= 1

            # M_t at t calculated from previous time step Mt and transitions thru matrix multiplication
            # .T is transpose
            # @ is matrix multiplication
            M_t[:,t] = transition_matrix_t.T @ Mt


        # Total productivity(?)
        Y_t                 = np.sum(M_t[[0,2,4,6,8,10,14]], axis=0) + \
                                self.A_rel * np.sum(M_t[[1,3,5,7,9,11,15]], axis=0)
        Reported_T_start    = self.pop * (tau_t + self.δ) * (M_t[4] + M_t[5])
        Reported_T_start[0] = 0
        Reported_T          = np.cumsum(Reported_T_start)

        Reported_D      = Reported_T[13::14] # Note: 13::14 refers to time indices, i.e. 'end of day for all days'
        Notinfected_D   = np.sum(M_t[[0,1,2,3]], axis=0)[13::14]
        Unreported_D    = np.sum(M_t[[4,5]], axis=0)[13::14]
        Infected_D      = np.sum(M_t[[12,13]], axis=0)[13::14]
        Infected_in_Q   = np.sum(M_t[[5,7,11,13]], axis=0)[13::14]      # includes all infected in quarantine including false negs
        Infected_not_Q  = np.sum(M_t[[4,6,10,12]], axis=0)[13::14]       # includes false negatives
        False_pos       = np.sum(M_t[[8,9]], axis=0)[13::14]
        False_neg       = np.sum(M_t[[10,11]], axis=0)[13::14]
        Recovered_D     = np.sum(M_t[[14,15]], axis=0)[13::14]
        Dead_D          = M_t[16][13::14]
        Infected_T      = np.sum(M_t[4:8], axis=0) + np.sum(M_t[10:16], axis=0)
        Y_D             = Y_t[13::14]

        return Reported_D, Notinfected_D, Unreported_D, Infected_D, \
                False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t


    def solve_model(self):
        Reported_D_base, Notinfected_D_base, Unreported_D_base, Infected_D_base, \
                False_pos_base, False_neg_base, Recovered_D_base, Dead_D_base, Infected_T_base, Infected_not_Q_base, Infected_in_Q_base, Y_D_base, M_t_base = \
                self.solve_case(self.baseline)
        Tstar = np.argwhere(Reported_D_base>100)[0][0]
        YearsPlot = 3
        Tplot = np.arange(Tstar, min(Tstar + YearsPlot * 365, self.T/self.Δ_time) + .5, 1)
        Xplot = np.arange(0, len(Tplot))
        self.Tstar = Tstar

        self.common_quarantine['d_start_exp'] = (Tstar+1) * self.Δ_time + \
                self.policy_offset * self.Δ_time

        Reported_D_com, Notinfected_D_com, Unreported_D_com, Infected_D_com, \
                False_pos_com, False_neg_com, Recovered_D_com, Dead_D_com, Infected_T_com, Infected_not_Q_com, Infected_in_Q_com, Y_D_com, M_t_com = \
                self.solve_case(self.common_quarantine)

        return Reported_D_com, Infected_D_com, Dead_D_com, Y_D_com, False_pos_com, False_neg_com, Infected_not_Q_com, Infected_in_Q_com

    def run_experiment(self, τ, Δ, test_sens, test_spec):

        τ_A_daily_target = τ

        r_U_daily_target	= 0
        r_N_daily_target	= 0
        r_P_daily_target	= 0
        r_R_daily_target	= self.r_high
        r_AP_daily_target   = self.r_AP # self.r

        ξ_U_daily_target   = self.ξ_base
        ξ_P_daily_target   = self.ξ_base_high
        ξ_N_daily_target   = self.ξ_base*Δ
        ξ_R_daily_target   = 0
        test_sens_exp      = test_sens
        test_spec_exp      = test_spec

        self.test_and_quarantine = {
            'τA'            : (1+τ_A_daily_target)**(1./self.Δ_time)-1,
            'test_sens'     : test_sens_exp,
            'test_spec'     : test_spec_exp,
            'ξ_U'           : (1+ξ_U_daily_target)**(1./self.Δ_time)-1,
            'ξ_P'           : (1+ξ_P_daily_target)**(1./self.Δ_time)-1,
            'ξ_N'           : (1+ξ_N_daily_target)**(1./self.Δ_time)-1,
            'ξ_R'           : (1+ξ_R_daily_target)**(1./self.Δ_time)-1,
            'r_U'           : (1+r_U_daily_target)**(1./self.Δ_time)-1,
            'r_P'           : (1+r_P_daily_target)**(1./self.Δ_time)-1,
            'r_AP'          : (1 + r_AP_daily_target) ** (1. / self.Δ_time) - 1,
            'r_N'           : (1+r_N_daily_target)**(1./self.Δ_time)-1,
            'r_R'           : (1+r_R_daily_target)**(1./self.Δ_time)-1,
            'experiment'    : "baseline_vaccine_tag"
        }

        self.test_and_quarantine['d_start_exp'] = (self.Tstar+1) * self.Δ_time + \
                self.policy_offset * self.Δ_time

        Reported_D_test, Notinfected_D_test, Unreported_D_test, Infected_D_test, \
                False_pos_test, False_neg_test, Recovered_D_test, Dead_D_test, Infected_T_test,  Infected_not_Q_test, Infected_in_Q_test, Y_D_test, M_t_test = \
                self.solve_case(self.test_and_quarantine)

        return Reported_D_test, Infected_D_test, Dead_D_test, Y_D_test, False_pos_test, False_neg_test, Infected_not_Q_test, Infected_in_Q_test,


def generate_plots(Δ, τ, test_sens, test_spec, ξ_base, A_rel, r_AP, d_vaccine, rel_ρ, δ_param, \
             ωR_param, π_D, R_0, rel_λ, initial_infect, slide_var):

    
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

    fig = make_subplots(4, 2, print_grid = False, \
                        subplot_titles=("A. Reported cases", "B. Current symptomatic cases", "C. Deaths - Cumulative", "D. Current output", "E. False positives", "F: False negatives", "G: Infected, not quarantined", "H: Infected, in quarantine"),
                        vertical_spacing = .2)

    print("Creating a corona model with sensitivity = ", test_sens, " and specificity = ", test_spec)
    model = corona_model(ξ_base, A_rel, r_AP, d_vaccine, rel_ρ, δ_param, \
                 ωR_param, π_D, R_0, rel_λ, initial_infect)

    Reported_D_com, Infected_D_com, Dead_D_com, Y_D_com, False_pos_com, False_neg_com, Infected_not_Q_com, Infected_in_Q_com, = model.solve_model()

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

    fig.add_scatter(y = Reported_D_com, row = 1, col = 1, visible = True, showlegend = True,
                    name = 'Common Quarantine', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Infected_D_com, row = 1, col = 2, visible = True, showlegend = False,
                    name = 'Common Quarantine', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Dead_D_com, row = 2, col = 1, visible = True, showlegend = False,
                    name = 'Common Quarantine', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y = Y_D_com, row = 2, col = 2, visible = True, showlegend = False,
                    name = 'Common Quarantine', line = dict(color = (colors[0]), width = 3, dash = styles[0]))
    fig.add_scatter(y=False_pos_com, row=3, col=1, visible=True, showlegend=False,
                    name='Common Quarantine', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=False_neg_com, row=3, col=2, visible=True, showlegend=False,
                    name='Common Quarantine', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=Infected_not_Q_com, row=4, col=1, visible=True, showlegend=False,
                    name='Common Quarantine', line=dict(color=(colors[0]), width=3, dash=styles[0]))
    fig.add_scatter(y=Infected_in_Q_com, row=4, col=2, visible=True, showlegend=False,
                    name='Common Quarantine', line=dict(color=(colors[0]), width=3, dash=styles[0]))

    if slide_var == 1: #Slide over τ
        prd = product(τ, [Δ], [test_sens], [test_spec])
        slider_vars = τ
        slider_varname = "τ"

    if slide_var == 2: #Slide over Δ
        prd = product([τ], Δ, [test_sens], [test_spec])
        slider_vars = Δ
        slider_varname = "Δ"

    if slide_var == 3:  # Slide over test_sens
        prd = product([τ], [Δ], test_sens, [test_spec])
        slider_vars = test_sens
        slider_varname = "test sensitivity"

    if slide_var == 4:  # Slide over test_spec
        prd = product([τ], [Δ], [test_sens], test_spec)
        slider_vars = test_spec
        slider_varname = "test specificity"

    pool = Pool(os.cpu_count())
    results = pool.starmap(model.run_experiment, prd)

    for j in range(len(slider_vars)):

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

    steps = []
    for i in range(len(slider_vars)):
        step = dict(
            method = 'restyle',
            args = [{'visible': ['legendonly'] * len(fig.data)},
                    {'showlegend': ['False'] * len(fig.data)}],
            label = slider_varname + ' = \n'+'{}'.format(round(slider_vars[i], 3))
        )
        step['args'][1]['showlegend'][0] = True
        step['args'][1]['showlegend'][8 + i * 8] = True
        for j in range(8):
            step['args'][0]['visible'][int(j)] = True
        for j in range(8):
            step['args'][0]['visible'][8 + j + i * 8] = True
        steps.append(step)

    sliders = [dict(
        steps = steps
    )]

    fig.layout.sliders = sliders
    for i in fig['layout']['annotations']:
        i['font'] = dict(color='black', size = 16)
    fig['layout'].update(height=800, width=1000, showlegend = False)

    fig['layout']['xaxis1'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 600], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis2'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 600], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis3'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 600], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis4'].update(title = go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range = [0, 600], \
                                   gridcolor = 'rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis5'].update(title=go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 600], \
                                    gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis6'].update(title=go.layout.xaxis.Title(
                                text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 600], \
                                    gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis7'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 600], \
        gridcolor='rgb(220,220,220)', showline=True, linewidth=1, linecolor='black', mirror=True)
    fig['layout']['xaxis8'].update(title=go.layout.xaxis.Title(
        text='Days since 100th case (3/4/2020)', font=dict(color='black')), range=[0, 600], \
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

    # fig['layout']['margin'].update(l=70, r=70, t=20, b=70)

    fig['layout']['plot_bgcolor'] = 'rgba(0,0,0,0)'

    return fig

def generate_plots_2d(Δ, τ, test_sens, test_spec, ξ_base, A_rel, d_vaccine, rel_ρ, δ_param, \
             ωR_param, π_D, R_0, rel_λ, initial_infect):

    model = corona_model(ξ_base, A_rel, d_vaccine, rel_ρ, δ_param, \
                 ωR_param, π_D, R_0, rel_λ, initial_infect)

    Reported_D, Infected_D, Dead_D, Y_D, Reported_D_com, Infected_D_com, \
        Dead_D_com, Y_D_com = model.solve_model()

    prd = product(τ, Δ)

    pool = Pool(os.cpu_count())
    results = pool.starmap(model.run_experiment, prd)

    return Reported_D, Infected_D, Dead_D, Y_D, Reported_D_com, Infected_D_com, \
        Dead_D_com, Y_D_com, results, prd
