#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=150M
#SBATCH --array=0-4

module restore covid_opt

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='base_case_lockdown_opt_with_limited_sens075_general_testing' ;;
    1)  RUN='base_case_lockdown_opt_with_limited_spec075_general_testing' ;;
    2)  RUN='base_case_lockdown_opt_with_limited_sens090_general_testing' ;;
    3)  RUN='base_case_lockdown_opt_with_limited_spec090_general_testing' ;;
    4)  RUN='romer_R0_1.25'  ;;
    

esac

srun python Covid_Run_Optimizer.py $RUN
