# To Do
FOCUS: results
## Paper: 
- (Check comments)
- Check methodology section
- Update figures with new clustering
- Improve results and make sure conclusions are in sync
- Remove figure 6

- Figure 5 - which bits to select
	- no need to show false positives if testing included
	- top 4 figures with 3 representative (red, blue, green)
	- leave testing profiles out
- Add references e.g. from epidemic modeling (see folders)

- Add the ‘No control’ as a reference case.

- Vocabulary:
	+ “Reduction” vs. “Pruning”
		+ -> pruning... check their meaning?
	+ How much do we talk about COVID-19? Don’t be afraid.
	+ “Robustness” -> use and explain until better one found. Maybe there’s a better word? - “Distributional robustness” / “Probabilistic risk constraint”
	- How to refer to time element in the decision variable vector? (‘time series’ - NO, ‘time profile’....)?
	- Change “Solution space” to “Decision variable space”?
- “Main results”
	- Implementability
	- Robustness
	- Policy imiplications for COVID
		- couunterintuitivee imiperfect test -> less deaths
## Analysis:
+ Reorder analysis: risk elimination first (before non-domination check)
- WHAT WAS THE OPTIMIZATION SENSITIVITY CHECK...?
	- Perhaps e.g.95% specificity
- Try syncing result files thru github...?
- Explore 'winning policy optioins' and how the epidemic behaves
- Add sens 85% (spec 100%), specificity 85% (sens 100%)
- The number of clusters for c-r-c pruning.
- sample size 10k (for sensitivity)
+ Model different ICU capacities. Can be used e.g. 'make sure less than 10% overload P for best output solutions etc'
+ Cluster -> combine -> risk (check if front shape enough in line with combine before cluster)
+ Risk analysis for 3 most significant ones, separatetly and together... Separately for cognitive support for reader
+ Why no change in overload P as f of ICU capa.?
- (Restudy effects of T_rec (recovery time, i.e. time for terminal values))
+ Scale output to range [0,1]
	+ Check that results OK?
+ Run with imperfect s&s for T&T
- Compare romer sens spec 085, 095, 099 policies - why 0.95 has worst death performance?
- Start policies to include (full, full/2, none)?

- Long term: Strategic uncertainty
- Long term: trigger controls by epidemic state
## General:


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

Replace one solution with perfect testing alternative OR fix the clustering algo so that it keeps at least one from each different group.
Focus on facilitiatioin and arriviing at a decision
Can work with any model + optimization method
Requires little interaction, can be quite fastly implemented (process implementability)
Robert J. Lempert (RDM)
Include schematic of the process
Include summary figure: obj + risk + solutions

Prepare runs to answer what is the effect of testing capacity constraints on risk and results
Update objectives (Excess hospital capacity & econ., keep deaths available.)
Just deaths + output as objectives.
Include optimization side in the default params
Check terminal cost vs. aggergated cost - see if reasonable
Risk measure: use 'probability to exceed capacity at teach time step' -> time vectors of binary overload incidcators -> combine samples to probabilities. Maximum P below treshold, e.g. 5% or lower. ('It is unacceptable for the probability to be higher than X at any time)
Uncertainty about initial infected (to model effects of delays).
Triton calculations for risk analysis

Solidify result collection, especially ensure order of solutions and obj values is maintained and can be sorted based on criteria (NOTE: pandas maintains order -> save everything to 1 dataframe per run -> add risk levels and sort according to risk level -> exclude too risk or e.g. 40 most risky ones (study) -> cluster.
Study and compare risk analysis and clustering output. Risky ones in all clusters or only in few?
Result presentation
Add base lambda as an argument for runs.
Include all main sensitivity parameters to risk analysis
Implement constraints for risk analysis.
Try risk analysis first with hospital capacity constraints / ranking - per TIME STEP - binary (exceeded / not exceeded)
Add terminal costs for all objectives
Initial conditions
Parameter value literature review
Only output as objective, costs recorded. -> 
Terminal cost sensitivity analysis
Fix testing constraints

- Constrain number of test through a NSGA-II constraint instead of test rate, as test rate isn't a 'stable' measure of number of tests.? Parallel risk model optimization
Salvage value, e.g. cost for number of infected.
Check time horizon sensitivity
Tradeoff between testing quality and capacity
What drives infections?
Remove all known infected not quarantined compartments as unnecessary
Quarantine traced persons first, then test if possible?

Finalize epidemic model
Known uninfected to unknown uninfected 
(see Matthias memo on rate) sigma=1/7)
Test & trace (Lauri & Matthias)
Check British 'only testing' case: sens, spec, rate, costs
Think through most important sensitivities:
Initial conditions
Parameter values: R_0, test_sensitivities, testing resources, 
Second shock
Delays: when implementation starts
Parallel runs on Triton
Adaptive Monte Carlo for sensitivity combinations
Make sure clustering works well in combo cases. This probably requires scaling the different control type values to similar scale.
Switch betat distr for pii_D Beta(1.45,95)
Control times to 1/month (30 days)
Lockdown to 0.2
Population to 100M
N_ICU scale according tot population
Max daily tests: divide by 3
sigma_param sensitivity (bit lower?)
Quarantine lambda to lower values: 0.1 (lower bound from Berger)
set testing ratet UL high enough
Sensitivities for N_ICU as separate scenarios

Powerpoint slides of results (when ready)
Simulate some worst case samples and see how they can be that bad!
Parameter values and extreme results
Soft lockdown with test and trace
Role of base lambda in risk results

COVID <-> dynamic, complex, how is it hard to make a decision?
Sharpen the story ono implementability and decision space considerations
Vaccination - mention somehow
Contributions
	 make more generic: procedure
	 applied side: 
	 Currrent ‘Thirdly’: what the effects of diffrent controls in MOO context..
	 Matthias suggested: 1st process. ‘sense making’
	 Robustness: highlight more (even own paragraph), note MORDM
Complexity of the model (medium)
	 Many realistic models are non-linear...
	 Model itself
	inclusion of sensitivity analysis, as an example, we’ll focus on COVID
EJOR style files
Insight: some parts of the Pareto front contain solutions where epidemic picks up again
Insight: knee solutions: characteristics
Insight: imperfect testing substituets lockdown
Discussion: smarter sampling techniques
Result cases: LD, testing, LD+testing, LD+T&T
Include Multiobjetive optimization as 1st contribution with also postprocessing (2nd) and 'insights' as 3rd.

Work with Combine - Risk - Cluster approach
# Maybe

Add 'perfect_testing_basecase' for background. 

Time step sensitivity analysis: how many control times and with what intervals?
