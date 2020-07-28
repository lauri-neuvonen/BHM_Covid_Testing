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

from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter


from policy_epidemic_model_code import *


termination = get_termination("n_gen", 200)

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
                         xl=np.array([0,0,0,0]),
                         xu=np.array([0.1,0.1,0.1,0.1]))

    def _evaluate(self, x, out, *args, **kwargs):

        f1 = []
        f2 = []
        for j in range(len(x[:,1])):  # evaluate f1 and f2 for each individual
            # TODO:  implement calculation of f1, f2 for all individuals (represented by rows in x)

            policy = {} # this will hold the policy in format suitable for input to the epidemic model
            for (i, t) in enumerate(self.policy_control_days):
                policy[t] = x[j,i]

            Reported_D, Notinfected_D, Unreported_D, Infected_D, \
            False_pos, False_neg, Recovered_D, Dead_D, Infected_T, Infected_not_Q, Infected_in_Q, Y_D, M_t, Y_total \
                = self.model.solve_case(self.model.baseline, policy)

            # objectives scaled to roughly same scale
            f1.append(Dead_D[-1]*1000000)
            f2.append(-Y_total/10000) # algorithm minimizes, Y_total needs to be max'd -> negative

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

problem = COVID_policy(corona_model, policy_control_days)

algorithm = NSGA2(
    pop_size=10,
    n_offsprings=6,
    sampling=get_sampling("real_random"),
    crossover=get_crossover("real_sbx", prob=0.9, eta=15),
    mutation=get_mutation("real_pm", eta=20),
    eliminate_duplicates=True
)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               #pf=problem.pareto_front(use_cache=False),
               save_history=True,
               verbose=True)

with open('optimization_result.txt', 'w') as f:
    print(res.X, file=f)

with open('optimization_objectives.txt', 'w') as f:
    print(res.F, file=f)


# get the pareto-set and pareto-front for plotting
#ps = problem.pareto_set(use_cache=False, flatten=False)
#pf = problem.pareto_front(use_cache=False, flatten=False)

# Design Space
plot = Scatter(title = "Design Space", axis_labels="x")
plot.add(res.X, s=30, facecolors='none', edgecolors='r')
#plot.add(ps, plot_type="line", color="black", alpha=0.7)
plot.do()
#plot.apply(lambda ax: ax.set_xlim(0.0, 5.0))
#plot.apply(lambda ax: ax.set_ylim(0.0, 0.2))
plot.show()

# Objective Space
plot = Scatter(title = "Objective Space")
plot.add(res.F)
plot.do()
#plot.add(pf, plot_type="line", color="black", alpha=0.7)
plot.show()


import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume

# create the performance indicator object with reference point (4,4)
metric = Hypervolume(ref_point=np.array([1.0, 1.0]))

# collect the population in each generation
pop_each_gen = [a.pop for a in res.history]

# receive the population in each generation
obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]

# calculate for each generation the HV metric
hv = [metric.calc(f) for f in obj_and_feasible_each_gen]

# function evaluations at each snapshot
n_evals = np.array([a.evaluator.n_eval for a in res.history])

# visualze the convergence curve
plt.plot(n_evals, hv, '-o')
plt.title("Convergence")
plt.xlabel("Function Evaluations")
plt.ylabel("Hypervolume")
plt.ylim([0.0,1.0])
plt.show()