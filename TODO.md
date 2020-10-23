# To Do

## General:
- Add terminal costs for all objectives
- Only output as objective, costs recorded. -> 
- Terminal cost sensitivity analysis
- Fix testing constraints
- 
- - Constrain number of test through a NSGA-II constraint instead of test rate, as test rate isn't a 'stable' measure of number of tests.? Parallel risk model optimization
- Salvage value, e.g. cost for number of infected.
- Check time horizon sensitivity
- Tradeoff between testing quality and capacity
- What drives infections?
+ Remove all known infected not quarantined compartments as unnecessary
+ Quarantine traced persons first, then test if possible?

- Finalize epidemic model
+ Known uninfected to unknown uninfected <br>(see Matthias memo on rate) sigma=1/7)
+ Test & trace (Lauri & Matthias)
- Check British 'only testing' case: sens, spec, rate, costs
+ Think through most important sensitivities:
* Initial conditions
* Parameter values: R_0, test_sensitivities, testing resources, 
- Second shock
+ Delays: when implementation starts
+ Parallel runs on Triton
- Adaptive Monte Carlo for sensitivity combinations

## Philosophical
- Difference measures for Pareto-optimal solutions...
## Epidemic modeling:

+  Create cost objective
- Scale lambda for output? Research literature for lockdown strength -> output change
* Think through the epidemic model and see what improvements could be made.
+ Recovery without symptoms - done thru retesting for FP
- Think about false positives: Now they mostly stay as false positives until vaccination as r^+ = 0. Some of them get sick and exit the state like that. Should they maybe have a release rate back to NA, NQ? Because: If you're identified as false positive but don't develop symptoms for a long time, you'd probably be retested or released as 'healed' (even though you'd not have immunity). This large amount of quarantined false positives creates a large reduction effect on the epidemic.

- Retesting for Paul Romer case: Change transitions from known uninfected's -> unknown uninfected + 

# Done

Implement capability to set imperfect testing params through notebook
1) Study why testing makes Infected Symptomatic compartments grow earlier - A: errors and bugs

Implement sliders for sensitivity and specificity

# Maybe

Add 'perfect_testing_basecase' for background. 
