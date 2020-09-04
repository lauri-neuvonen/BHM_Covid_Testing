#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-3

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='base_case_lockdown_opt' ;;
    1)  RUN='base_case_lockdown_opt_R0_4.0' ;;

esac

srun python Covid_Run_Optimizer.py $RUN