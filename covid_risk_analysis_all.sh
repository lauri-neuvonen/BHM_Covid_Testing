#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=110M
#SBATCH --array=0-13

module restore covid_opt

SAMPLE_SIZE=1000

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='romer' ;;
    1)  RUN='romer_spec_085' ;;
    2)  RUN='romer_sens_085' ;;
    3)  RUN='romer_sens_spec_085' ;;
    4)  RUN='romer_sens_spec_095' ;;
    5)  RUN='romer_sens_spec_099' ;;
    6)  RUN='base_case_lockdown_opt' ;;
    7)  RUN='base_case_lockdown_opt_with_limited_general_testing' ;;
    8)  RUN='base_case_lockdown_opt_with_limited_sens075_general_testing'  ;;
    9)  RUN='base_case_lockdown_opt_with_limited_spec075_general_testing'  ;;
    10)  RUN='base_case_lockdown_opt_with_limited_sens090_general_testing'  ;;
    11)  RUN='base_case_lockdown_opt_with_limited_spec090_general_testing'  ;;
    12)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing' ;;
    13)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing'  ;;

esac

srun python risk_analysis.py $SAMPLE_SIZE $RUN
