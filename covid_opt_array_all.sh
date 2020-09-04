#!/bin/bash
#SBATCH --time=0:24:00
#SBATCH --mem-per-cpu=200M
#SBATCH --array=0-30

module restore covid_opt

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='romer' ;;
    1)  RUN='romer_R0_4.0' ;;
    2)  RUN='romer_R0_4.0_sens_spec_075' ;;
    3)  RUN='romer_6d_incubation' ;;
    4)  RUN='romer_8d_incubation' ;;
    5)  RUN='romer_8d_incubation_sens_spec_075' ;;
    6)  RUN='romer_3d_delay' ;;
    7)  RUN='romer_7d_delay' ;;
    8)  RUN='romer_14d_delay' ;;
    9)  RUN='romer_sens_spec_075' ;;
    10)  RUN='romer_spec_075' ;;
    11)  RUN='romer_sens_075' ;;
    12)  RUN='romer_sens_spec_085' ;;
    13)  RUN='romer_sens_spec_090' ;;
    14)  RUN='romer_sens_spec_095' ;;
    15)  RUN='romer_R0_4.0_sens_spec_085' ;;
    16)  RUN='romer_R0_4.0_sens_spec_090' ;;
    17)  RUN='romer_R0_4.0_sens_spec_095' ;;
    18)  RUN='base_case_lockdown_opt' ;;
    19)  RUN='base_case_lockdown_opt_R0_4.0' ;;
    20)  RUN='base_case_lockdown_opt_with_limited_general_testing' ;;
    21)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing' ;;
    22)  RUN='base_case_lockdown_opt_3d_delay' ;;
    23)  RUN='base_case_lockdown_opt_7d_delay' ;;
    24)  RUN='base_case_lockdown_opt_14d_delay' ;;
    25)  RUN='base_case_lockdown_opt_28d_delay' ;;
    26)  RUN='base_case_6d_incubation' ;;
    27)  RUN='base_case_8d_incubation' ;;
    28)  RUN='base_case_lockdown_opt' ;;
    29)  RUN='base_case_lockdown_opt_R0_4.0' ;;
    30)  RUN='romer_28d_delay'  ;;


esac

srun python Covid_Run_Optimizer.py $RUN
