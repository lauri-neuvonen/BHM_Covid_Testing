import policy_epidemic_model_code
import numpy as np

from run_tools import create_epidemic_model, Policy, Policy_template
from run_definitions import get_runs_definitions
import argparse

parser = argparse.ArgumentParser(description='Run optimization run using an NSGA-II algorithm for COVID-19 simulator')
parser.add_argument('sample_size', type=int, help='sample size for risk analysis')
parser.add_argument('runs', type=str, nargs='+', help='dictionary of run values to be changed from defaults')

args = parser.parse_args()
sample_size = args.sample_size # sample size from input
run_list = args.runs # runs to be analysed

run_definitions = get_runs_definitions()

runs = { run: run_definitions[run] for run in run_list } # filters the correct run definitions based on given run list

for run in runs:
    print("run: ", runs[run])
    sample_list = np.zeros(sample_size)

    analysis_params = {}

    analysis_params['R_0']= {
        '_dist': np.random.normal,
        'exp': runs[run]['R_0'],
        'var': 0.2,
    }

    analysis_params['pii_D'] = {
        '_dist': np.random.normal,
        'exp': runs[run]['pii_D'],
        'var': 0.02
    }

    for i in range(0, sample_size):
        sample_instance = {}
        for param in analysis_params:
            sample_instance[param] = analysis_params[param]['_dist'](analysis_params['exp'], analysis_params['var'])

   #debug
    print("samples: ", sample_list)






