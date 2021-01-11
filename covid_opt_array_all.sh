#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --mem-per-cpu=110M
#SBATCH --array=0-63

module restore covid_opt

MAX_GEN=3000

case $SLURM_ARRAY_TASK_ID in

    0)  RUN='romer' ;;
    1)  RUN='romer_R0_1.25'  ;;
    2)  RUN='romer_R0_4.0' ;;
    3)  RUN='romer_R0_4.0_sens_spec_075' ;;
    4)  RUN='romer_R0_4.0_sens_spec_085' ;;
    5)  RUN='romer_R0_4.0_sens_spec_090' ;;
    6)  RUN='romer_R0_4.0_sens_spec_095' ;;
    7)  RUN='romer_6d_incubation' ;;
    8)  RUN='romer_8d_incubation' ;;
    9)  RUN='romer_8d_incubation_sens_spec_075' ;;
    10)  RUN='romer_8d_incubation_sens_spec_090'  ;;
    11)  RUN='romer_3d_delay' ;;
    12)  RUN='romer_28d_delay'  ;;
    13)  RUN='romer_spec_085' ;;
    14)  RUN='romer_sens_085' ;;
    15)  RUN='romer_sens_spec_075' ;;
    16)  RUN='romer_sens_spec_085' ;;
    17)  RUN='romer_sens_spec_090' ;;
    18)  RUN='romer_sens_spec_095' ;;
    19)  RUN='romer_sens_spec_098' ;;
    20)  RUN='romer_sens_spec_099' ;;
    21)  RUN='romer_terminal_0.25' ;;
    22)  RUN='romer_terminal_1.0' ;;
    23)  RUN='base_case_lockdown_opt' ;;
    24)  RUN='base_case_lockdown_opt_terminal_0.25' ;;
    25)  RUN='base_case_lockdown_opt_terminal_1.0' ;;
    26)  RUN='base_case_lockdown_opt_R0_4.0' ;;
    27)  RUN='base_case_lockdown_opt_with_limited_general_testing' ;;
    28)  RUN='base_case_lockdown_opt_with_limited_sens075_general_testing'  ;;
    29)  RUN='base_case_lockdown_opt_with_limited_spec075_general_testing'  ;;
    30)  RUN='base_case_lockdown_opt_with_limited_sens090_general_testing'  ;;
    31)  RUN='base_case_lockdown_opt_with_limited_spec090_general_testing'  ;;
    32)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.75)_general_testing' ;;
    33)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.85)_general_testing'  ;;
    34)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.90)_general_testing'  ;;
    35)  RUN='base_case_lockdown_opt_with_limited_imperfect(0.95)_general_testing'  ;;
    36)  RUN='base_case_lockdown_opt_3d_delay' ;;
    37)  RUN='base_case_lockdown_opt_7d_delay' ;;
    38)  RUN='base_case_lockdown_opt_14d_delay' ;;
    39)  RUN='base_case_lockdown_opt_28d_delay' ;;
    40)  RUN='base_case_6d_incubation' ;;
    41)  RUN='base_case_8d_incubation' ;;
    42)  RUN='test_and_trace_lockdown_opt_eta50'  ;;
    43)  RUN='test_and_trace_lockdown_opt_eta75'  ;;
    44)  RUN='test_and_trace_lockdown_opt_eta95'  ;;
    45)  RUN='test_and_trace_lockdown_opt_eta100'  ;;
    46)  RUN='test_and_trace_lockdown_opt_eta50_R04'  ;;
    47)  RUN='test_and_trace_lockdown_opt_eta75_R04'  ;;
    48)  RUN='test_and_trace_lockdown_opt_eta100_R04' ;;
    49)  RUN='test_and_trace_lockdown_opt_eta50_R04_delta10'  ;;
    50)  RUN='test_and_trace_lockdown_opt_eta75_R04_delta10'  ;;
    51)  RUN='test_and_trace_lockdown_opt_eta100_R04_delta10'  ;;
    52)  RUN='combo_base_case' ;;
    53)  RUN='combo_base_case_R0_4.0'  ;;
    54)  RUN='combo_base_case_test_and_trace'  ;;
    55)  RUN='combo_sens_spec_0.95' ;;
    56)  RUN='combo_sens_spec_0.85' ;;
    57)  RUN='combo_R0_4.0_sens_spec_0.95'  ;;
    58)  RUN='combo_R0_4.0_sens_spec_0.85'  ;;
    59)  RUN='romer_no_limit'  ;;
    60)  RUN='romer_sens_spec_085_hi_mut_param' ;;
    61)  RUN='romer_sens_spec_085_low_mut_param' ;;
    62)  RUN='romer_sens_spec_085_hi_cross_param' ;;
    63)  RUN='romer_sens_spec_085_low_cross_param' ;;


esac

srun python Covid_Run_Optimizer.py $MAX_GEN $RUN