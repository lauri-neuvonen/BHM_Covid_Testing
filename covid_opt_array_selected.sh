#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --mem-per-cpu=110M
#SBATCH --array=0-58

module restore covid_opt

MAX_GEN=5000
SUFFIX='_f-cand'

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='combo_base_case' ;;
    1)  RUN='combo_base_case_tc_50000000_sens_spec_085' ;;
    2)  RUN='romer' ;;
    3)  RUN='base_case_lockdown_opt' ;;
    4)  RUN='test_and_trace_lockdown_opt_eta100' ;;
    5)  RUN='romer_tc_50000000_sens_spec_085' ;;

esac

srun python Covid_Run_Optimizer.py $MAX_GEN $RUN --file_suffix=$SUFFIX
