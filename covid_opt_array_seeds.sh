#!/bin/bash
#SBATCH --time=00:24:00
#SBATCH --mem-per-cpu=200M
#SBATCH --array=0-1

module restore covid_opt

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='base_case_lockdown_opt' ;;
    1)  RUN='romer' ;;

esac

srun python Covid_Run_Optimizer.py $RUN
