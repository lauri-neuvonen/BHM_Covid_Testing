#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=110M
#SBATCH --array=0-5

module restore covid_opt

SAMPLE_SIZE=1000
SET='full_results_f-cand' # 'medoid' or 'full_results'
SUFFIX=''

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='combo_base_case' ;;
    1)  RUN='combo_base_case_tc_50000000_sens_spec_085' ;;
    2)  RUN='romer' ;;
    3)  RUN='base_case_lockdown_opt' ;;
    4)  RUN='test_and_trace_lockdown_opt_eta100' ;;
    5)  RUN='romer_tc_50000000_sens_spec_085' ;;
    
esac


srun python risk_analysis.py $SAMPLE_SIZE $SET $RUN 
