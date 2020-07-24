### This program optimizes an SEIR epidemic model using NSGA-II
### NSGA implementation through pymoo
### Epidemic model is an extended version of the BHM model (see README for details)

# License: MIT
# Copyright Lauri Neuvonen, July 2020

import numpy as np
from pymoo.util.misc import stack
from pymoo.model.problem import Problem

from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation

from pymoo.factory import get_termination

from policy_epidemic_model_code import *

termination = get_termination("n_gen", 40)

ξ_base = 1.0
A_rel = 0.5
r_AP = 0.0
d_vaccine = 500
rel_ρ = 1.0
δ_param = 6
ωR_param = 14
π_D =  0.01
R_0 = 2.0
rel_λ = 0.5
initial_infect = 300

corona_model = optimizable_corona_model(ξ_base, A_rel, r_AP, d_vaccine, rel_ρ, δ_param, \
                 ωR_param, π_D, R_0, rel_λ,initial_infect)

policy_control_days = [15, 30, 45, 100]

class COVID_policy(Problem):

    def __init__(self, model, policy_control_days):

        self.model = model
        self.policy_control_days = policy_control_days

        super().__init__(n_var=len(policy_control_days),
                         n_obj=2,
                         n_constr=0,
                         xl=np.array([0,0]),
                         xu=np.array([1,1]))

    def _evaluate(self, x, out, *args, **kwargs):

        for i in len(x[:,1]):
            # TODO:  implement calculation of f1, f2 for all individuals (represented by rows in x)
            policy_control_levels = x[i,:]


            policy = {}
            for (i, t) in enumerate(policy_control_days):
                policy[t] = x[:,i]

        Reported_D, Notinfected_D, Unreported_D, Infected_D, \
        False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total \
        = self.model.solve_policy_case(self.model.baseline, policy)
        ## For reference - TODO delete when not needed:
        #f1 = x[:,0]**2 + x[:,1]**2
        #f2 = (x[:,0]-1)**2 + x[:,1]**2

        #g1 = 2*(x[:, 0]-0.1) * (x[:, 0]-0.9) / 0.18
        #g2 = - 20*(x[:, 0]-0.4) * (x[:, 0]-0.6) / 4.8


        out["F"] = np.column_stack([f1, f2])
        #out["G"] = np.column_stack([g1, g2])


    # --------------------------------------------------
    # Pareto-front - not necessary but used for plotting
    # --------------------------------------------------
    def _calc_pareto_front(self, flatten=True, **kwargs):
        f1_a = np.linspace(0.1**2, 0.4**2, 100)
        f2_a = (np.sqrt(f1_a) - 1)**2

        f1_b = np.linspace(0.6**2, 0.9**2, 100)
        f2_b = (np.sqrt(f1_b) - 1)**2

        a, b = np.column_stack([f1_a, f2_a]), np.column_stack([f1_b, f2_b])
        return stack(a, b, flatten=flatten)

    # --------------------------------------------------
    # Pareto-set - not necessary but used for plotting
    # --------------------------------------------------
    def _calc_pareto_set(self, flatten=True, **kwargs):
        x1_a = np.linspace(0.1, 0.4, 50)
        x1_b = np.linspace(0.6, 0.9, 50)
        x2 = np.zeros(50)

        a, b = np.column_stack([x1_a, x2]), np.column_stack([x1_b, x2])
        return stack(a,b, flatten=flatten)

problem = MyProblem()

algorithm = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=get_sampling("real_random"),
    crossover=get_crossover("real_sbx", prob=0.9, eta=15),
    mutation=get_mutation("real_pm", eta=20),
    eliminate_duplicates=True
)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               pf=problem.pareto_front(use_cache=False),
               save_history=True,
               verbose=True)