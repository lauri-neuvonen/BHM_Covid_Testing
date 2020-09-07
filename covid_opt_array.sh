#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=150M
#SBATCH --array=0-4

module restore covid_opt

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='base_case_lockdown_opt_28d_delay' ;;
    1)  RUN='romer_28d_delay' ;;
    2)  RUN='romer_sens_075'  ;;
    3)  RUN='romer_spec_075'  ;;
    4)  RUN='romer_R0_1.25'  ;;
    

esac

srun python Covid_Run_Optimizer.py $RUN
