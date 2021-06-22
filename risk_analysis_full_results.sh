#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=110M
#SBATCH --array=0-80

module restore covid_opt

SAMPLE_SIZE=1000
SET='full_results_f-cand' # 'medoid' or 'full_results'
SUFFIX=''

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
    14)  RUN='combo_base_case' ;;
    15)  RUN='romer_tc_1000000_sens_spec_085' ;;
    16)  RUN='romer_tc_2500000_sens_spec_085' ;;
    17)  RUN='romer_tc_5000000_sens_spec_085' ;;
    18)  RUN='romer_tc_25000000_sens_spec_085' ;;
    19)  RUN='romer_tc_50000000_sens_spec_085' ;;
    20)  RUN='romer_tc_100000000_sens_spec_085' ;;
    21)  RUN='romer_tc_1000000' ;;
    22)  RUN='romer_tc_2500000' ;;
    23)  RUN='romer_tc_5000000' ;;
    24)  RUN='romer_tc_25000000' ;;
    25)  RUN='romer_tc_50000000' ;;
    26)  RUN='romer_tc_100000000' ;;
    27)  RUN='combo_base_case_test_and_trace' ;;
    28)  RUN='combo_base_case_tc_1000000' ;;
    29)  RUN='combo_base_case_tc_2500000' ;;
    30)  RUN='combo_base_case_tc_5000000' ;;
    31)  RUN='combo_base_case_tc_25000000' ;;
    32)  RUN='combo_base_case_tc_50000000' ;;
    33)  RUN='combo_base_case_tc_100000000' ;;
    34)  RUN='test_and_trace_lockdown_opt_eta50' ;;
    35)  RUN='test_and_trace_lockdown_opt_eta75' ;;
    36)  RUN='test_and_trace_lockdown_opt_eta100' ;;
    37)  RUN='combo_base_case_test_and_trace_tc1000000' ;;
    38)  RUN='combo_base_case_test_and_trace_tc2500000' ;;
    39)  RUN='combo_base_case_test_and_trace_tc5000000' ;;
    40)  RUN='combo_base_case_test_and_trace_tc25000000' ;;
    41)  RUN='combo_base_case_test_and_trace_tc50000000' ;;
    42)  RUN='combo_base_case_test_and_trace_tc100000000' ;;
    43)  RUN='combo_base_case_test_and_trace_ss085' ;;
    44)  RUN='combo_base_case_test_and_trace_tc1000000_ss085' ;;
    45)  RUN='combo_base_case_test_and_trace_tc2500000_ss085' ;;
    46)  RUN='combo_base_case_test_and_trace_tc5000000_ss085' ;;
    47)  RUN='combo_base_case_test_and_trace_tc25000000_ss085' ;;
    48)  RUN='combo_base_case_test_and_trace_tc50000000_ss085' ;;
    49)  RUN='combo_base_case_test_and_trace_tc100000000_ss085' ;;
    50)  RUN='combo_base_case_sens_spec_085' ;;
    51)  RUN='combo_base_case_tc_500000_sens_spec_085' ;;
    52)  RUN='combo_base_case_tc_1000000_sens_spec_085' ;;
    53)  RUN='combo_base_case_tc_2500000_sens_spec_085' ;;
    54)  RUN='combo_base_case_tc_5000000_sens_spec_085' ;;
    55)  RUN='combo_base_case_tc_25000000_sens_spec_085' ;;
    56)  RUN='combo_base_case_tc_50000000_sens_spec_085' ;;
    57)  RUN='test_and_trace_lockdown_opt_eta100_sens_spec_085' ;;
    58)  RUN='test_and_trace_lockdown_opt_eta50_sens_spec_085' ;;
    59)  RUN='romer --params R_0 --file_suffix=$_R_0_only' ;;
    60)  RUN='combo_base_case --params R_0 --file_suffix=$_R_0_only' ;;
    61)  RUN='base_case_lockdown_opt --params R_0 --file_suffix=$_R_0_only' ;;
    62)  RUN='romer_tc_50000000_sens_spec_085 --params R_0 --file_suffix=$_R_0_only' ;;
    63)  RUN='combo_base_case_tc_50000000_sens_spec_085 --params R_0 --file_suffix=$_R_0_only' ;;
    64)  RUN='test_and_trace_lockdown_opt_eta100 --params R_0 --file_suffix=$_R_0_only' ;;
    59)  RUN='romer --params delta_param --file_suffix=$_delta_only' ;;
    60)  RUN='combo_base_case --params delta_param --file_suffix=$_delta_only' ;;
    61)  RUN='base_case_lockdown_opt --params delta_param --file_suffix=$_delta_only' ;;
    62)  RUN='romer_tc_50000000_sens_spec_085 --params delta_param --file_suffix=$_delta_only' ;;
    63)  RUN='combo_base_case_tc_50000000_sens_spec_085 --params delta_param --file_suffix=$_delta_only' ;;
    64)  RUN='test_and_trace_lockdown_opt_eta100 --params delta_param --file_suffix=$_delta_only' ;;
    59)  RUN='romer --params pii_D --file_suffix=$_pii_D_only' ;;
    60)  RUN='combo_base_case --params pii_D --file_suffix=$_pii_D_only' ;;
    61)  RUN='base_case_lockdown_opt --params pii_D --file_suffix=$_pii_D_only' ;;
    62)  RUN='romer_tc_50000000_sens_spec_085 --params pii_D --file_suffix=$_pii_D_only' ;;
    63)  RUN='combo_base_case_tc_50000000_sens_spec_085 --params pii_D --file_suffix=$_pii_D_only' ;;
    64)  RUN='test_and_trace_lockdown_opt_eta100 --params pii_D --file_suffix=$_pii_D_only' ;;

esac


srun python risk_analysis.py $SAMPLE_SIZE $SET $RUN 
