import policy_epidemic_model_code
import numpy as np

from run_tools import create_epidemic_model, create_policy_control, Policy, Policy_template
from run_definitions import get_runs_definitions
import argparse
from run_definitions import (ksi_base_default, A_rel_default, r_AP_default, r_U_default, r_P_default, r_N_default, d_vaccine_default, rel_rho_default, delta_param_default,
    omegaR_param_default, pii_D_default, R_0_default, rel_lambda_param_default, gamma_param_default, initial_infect_default, daily_testing_rate_default,
    testing_sensitivity_default, testing_specificity_default, tau_TT_daily_default, eta_default, unknown_q_rate_default, recovered_q_rate_default,
    negative_q_rate_default, positive_q_rate_default, testing_cost_default)

parser = argparse.ArgumentParser(description='Run optimization run using an NSGA-II algorithm for COVID-19 simulator')
parser.add_argument('sample_size', type=int, help='sample size for risk analysis')
parser.add_argument('runs', type=str, nargs='+', help='dictionary of run values to be changed from defaults')

args = parser.parse_args()
sample_size = args.sample_size # sample size from input
run_list = args.runs # runs to be analysed

run_definitions = get_runs_definitions()

runs = { run: run_definitions[run] for run in run_list } # filters the correct run definitions based on given run list

### INPUT RISK ANALYSIS DEFINITIONS HERE! ###
### MAKE SURE CORRECT SOURCE FOR EXP is used ###

for run in runs:

    # get policies:

    run_policies_df = pd.read_csv('active_results/' + run + '_results.csv', delimiter=',')
    run_policies = run_policies_df.to_numpy()


    # get parameter value samples:
    print("run: ", runs[run])
    sample_list = []

    analysis_params = {}

    try:
        analysis_params['R_0']= {
            '_dist': np.random.lognormal,
            'dist_params': (np.log(runs[run]['R_0']**2/np.sqrt((runs[run]['R_0']**2 + 0.2))), np.log(1 + 0.2/runs[run]['R_0']**2))  # [mean, sigma] for lognormal
        }
    except KeyError:
        analysis_params['R_0'] = {
            '_dist': np.random.lognormal,
            'dist_params': (np.log(R_0_default**2/np.sqrt((R_0_default**2 + 0.2))), np.log(1 + 0.2/R_0_default**2))   # [mean, sigma] for lognormal
        }

    analysis_params['pii_D'] = {
        '_dist': np.random.lognormal,
        'dist_params': (
        np.log(pii_D_default ** 2 / np.sqrt((pii_D_default ** 2 + 0.0002))), np.log(1 + 0.0002 / pii_D_default ** 2))
        # [mean, sigma] for lognormal
    }

    for i in range(0, sample_size):
        sample_instance = {}
        for param in analysis_params:
            sample_instance[param] = analysis_params[param]['_dist'](*analysis_params[param]['dist_params'])
        sample_list.append(sample_instance)
   #debug
    print("samples: ", sample_list)

    for sample in sample_list:
        sample_run_params = runs[run].copy() # copies the original run (e.g. 'romer') for updating with sample values
        sample_run_params.update(sample)
        #debug:
        print("sample params: ", sample_run_params)

        model, model_case = create_epidemic_model(sample_run_params)
        policy = create_policy_control() # how was this done...





