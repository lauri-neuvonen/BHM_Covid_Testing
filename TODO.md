# To Do

## General:
### Decide on retesting / returning people to NANQ state / how transitions to infected states are modeled
In the original model it was assumed that once tested, perfect information on a persons state is always available.

This was a problem for imperfect testing, so sensitivity ans specificity corrections were introduced to these transitions.

If we return people to NANQ state from NANQ* (so information dissipates) this transition from NANQ* to IANQ* (and Q) becomes incoherent with the story. 

#### option 1: No transitions from known back to unknown states
#### option 2: Lose information
* Remove transitions from NANQ* and ...Q* to IANQ* and ...Q* - instead move to unknown infected states. Transitions from known uninfected states to known infected states then through unknown infected states.
* Also consider retesting Q* states to (test to release).

### other
- Known uninfected to unknown uninfected <br>(see Matthias memo on rate) sigma=1/7)
- Test & trace (Lauri & Matthias)
- Check British 'only testing' case: sens, spec, rate, costs
- Think through most important sensitivities:
* Initial conditions
* Parameter values: R_0, test_sensitivities, testing resources, 
- Second shock
- Delays: when implementation starts
+ Parallel runs on Triton
- Adaptive Monte Carlo for sensitivity combinations

## Philosophical
- Difference measures for Pareto-optimal solutions...
## Optimization:

- Add cost to objective visualization
- Change baseline in objective visualization to a line instead of bar
* Find good crossover and mutation parameters
* Create cost objective (first in epidemic model)

## Epidemic modeling:

+  Create cost objective
- Scale lambda for output? Research literature for lockdown strength -> output change
* Think through the epidemic model and see what improvements could be made.
- Recovery without symptoms
- Think about false positives: Now they mostly stay as false positives until vaccination as r^+ = 0. Some of them get sick and exit the state like that. Should they maybe have a release rate back to NA, NQ? Because: If you're identified as false positive but don't develop symptoms for a long time, you'd probably be retested or released as 'healed' (even though you'd not have immunity). This large amount of quarantined false positives creates a large reduction effect on the epidemic.

- Retesting for Paul Romer case: Change transitions from known uninfected's -> unknown uninfected + 

# Done

Implement capability to set imperfect testing params through notebook
1) Study why testing makes Infected Symptomatic compartments grow earlier - A: errors and bugs

Implement sliders for sensitivity and specificity

# Maybe

Add 'perfect_testing_basecase' for background. 
