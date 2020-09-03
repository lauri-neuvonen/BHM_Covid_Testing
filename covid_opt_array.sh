#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-3

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='romer' ;;
    1)  RUN='romer_R0_4.0'  ;;
    2)  RUN='romer_R0_4.0_sens_spec_075'  ;;
    3)  RUN='romer_6d_incubation'  ;;
esac

srun python Covid_Run_Optimizer.py $RUN