#!/bin/bash -i
# ==============================================================================

# Name: run_project.sh
# Author: Braden Yang
# Created: 01/09/23
# Description: Run all scripts in this project repository

# Notes
# - currently the simulation is only set to run on a single set of parameters and for
#   only 10 iterations
# - to run on more interations, it is advised to call multiple jobs, one for each
#   set of parameters in sim_param.txt

# ==============================================================================

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

usage() {
    echo "Usage: $0 [-i <int>] [-n <int>] [-h]"
    echo "Options:"
    echo "  -i <int>: only run step i"
    echo "  -n <int>: run all steps up to n (inclusive)"
    echo "  -h: Display this help message"
    echo ""
    echo "Details:"
    echo "  step 0: run 0_GetData"
    echo "  step 1: run 1_Harmonize"
    echo "  step 2: run 2_HeadToHead"
    echo "  step 3: run 3_Simulation (except visualizations)"
    echo "  step 4: run 3_Simulation (visualizations only)"
    echo "  step 5: run get_demographics.R, get_results_table.R, and get_excel.R"
    exit 1
}

# define defaults
run_up_to="1000"  # run everything by default
idx_to_run=""

# read arguments
while getopts i:n:h arg
do
	case $arg in
    i)  idx_to_run=${OPTARG};;
	n)	run_up_to=${OPTARG};;
    h)  usage;;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
	esac
done

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

stop_early () {
    if [[ ${run_up_to} -lt $1 ]]; then
        exit 0
    fi
}

run_step0 () {

    echo ""
    echo "++ Running 0_GetData ++"
    echo ""

    Rscript code/0_GetData/get_data.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data
    Rscript code/0_GetData/get_oasis_crossover.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data
    Rscript code/0_GetData/get_oasis_crossover.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data \
        --pvc

}

run_step1 () {

    echo ""
    echo "++ Running 1_Harmonize ++"
    echo ""

    Rscript code/1_Harmonize/centiloid.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data/centiloid
    python code/1_Harmonize/apply_combat.py \
        --data_dir ${proj_dir}/data \
        --odir ${proj_dir}/data/combat

    Rscript code/1_Harmonize/PEACE/apply_PEACE_nocovar.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data/head2head/PEACE
    Rscript code/1_Harmonize/PEACE/apply_PEACE_age_sex_apoe.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/data/head2head/PEACE

}

run_step2 () {

    echo ""
    echo "++ Running 2_HeadToHead ++"
    echo ""

    Rscript code/2_HeadToHead/crossover_analysis.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/figures/crossover
    Rscript code/2_HeadToHead/crossover_analysis.R \
        --wdir ${proj_dir} \
        --odir ${proj_dir}/figures/crossover \
        --plot_peace

}

run_step3 () {

    echo ""
    echo "++ Running 3_Simulation ++"
    echo ""

    # create amypos.csv to train ComBat on
    if [ ! -f ${sim_dir}/amypos.csv ]; then
        Rscript code/3_Simulation/1_simulate_data.R \
            -w ${proj_dir} \
            -o ${sim_dir} \
            -a
    fi

    # train ComBat
    if [ ! -d ${sim_dir}/trained_combat ]; then
        python code/3_Simulation/2_train_combat.py \
            -w ${proj_dir} \
            -i ${sim_dir}/amypos.csv \
            -o ${sim_dir}/trained_combat \
            -f ${feature_names_trunc}
    fi

    # train PEACE
    if [ ! -d ${sim_dir}/PEACE ]; then
        Rscript code/3_Simulation/PEACE/train_PEACE.R \
            -w ${proj_dir} \
            -o ${sim_dir}/PEACE \
            -m "nocovar"
        Rscript code/3_Simulation/PEACE/train_PEACE.R \
            -w ${proj_dir} \
            -o ${sim_dir}/PEACE \
            -m "age_sex_apoe"
    fi

    # run simulation loop
    i=1
    while [[ $i -le $n_param ]]; do

        param_arr=($(cat ${param_path} | awk -v line=${i} '{if (NR == line) print $0}'))
        n_subj="${param_arr[0]}"
        p_av45_placebo="${param_arr[1]}"
        p_av45_treatment="${param_arr[2]}"
        treatment_effect="${param_arr[3]}"
        sim_data_dir="${sim_dir}/${param_arr[4]}"

        echo "++ working on ${sim_data_dir} ++"

        # simulate data
        Rscript code/3_Simulation/1_simulate_data.R \
            -w ${proj_dir} \
            -o ${sim_data_dir} \
            -n ${n_simulation} \
            -p ${p_av45_placebo} \
            -q ${p_av45_treatment} \
            -s ${n_subj} \
            -t ${treatment_effect} \
            -r 42 \
            -x

        # harmonize simulated data
        # Centiloid, ComBat
        python code/3_Simulation/3_harmonize.py \
            -w ${proj_dir} \
            -i ${sim_data_dir}/sim.csv \
            -o ${sim_data_dir} \
            -c ${sim_dir}/trained_combat \
            -f ${feature_names_trunc}
        # PEACE, no covar
        Rscript code/3_Simulation/PEACE/apply_PEACE.R \
            -w ${proj_dir} \
            -o ${sim_data_dir} \
            -p ${sim_dir}/PEACE \
            -m "nocovar" \
            -f ${sim_data_dir}/sim.csv
        # PEACE, with covar
        Rscript code/3_Simulation/PEACE/apply_PEACE.R \
            -w ${proj_dir} \
            -o ${sim_data_dir} \
            -p ${sim_dir}/PEACE \
            -m "age_sex_apoe" \
            -f ${sim_data_dir}/sim.csv
        # longCombat
        Rscript code/3_Simulation/longComBat/apply_long_combat.R \
            -w ${proj_dir} \
            -o ${sim_data_dir} \
            -f ${sim_data_dir}/sim.csv
        
        # merge harmonization outputs
        Rscript code/3_Simulation/merge_harm.R \
            -w ${proj_dir} \
            -d ${sim_data_dir}

        # test for statistical significance
        for roi in $(cat ${feature_names_trunc}); do

            Rscript code/3_Simulation/4_statistical_test.R \
                -w ${proj_dir} \
                -i ${sim_data_dir}/sim_harm.csv \
                -o ${sim_data_dir} \
                -f ${roi} \
                -t

        done

        i=$(($i+1))
    
    done

    # collect all simulation results into one table
    if [ ! -f ${sim_dir}/stats.RDS ]; then
        Rscript code/3_Simulation/5_visualize.R \
            -w ${proj_dir} \
            -r \
            -d ${sim_dir} \
            -o ${sim_dir}
    fi

}

run_step4 () {

    # NOTE: this step uses the already-generated stats.RDS
    # to create the figures, since it takes a long time to
    # get stats.RDS from running 3_Simulation with 1000
    # iterations and multiple parameters (tracer mixing
    # proportion, treatment effect)

    echo ""
    echo "++ Running 3_Simulation (visualizations) ++"
    echo ""

    # create visualizations
    Rscript code/3_Simulation/5_visualize.R \
        -s ${sim_dir}/stats.RDS \
        -o figures/simulation
    Rscript code/3_Simulation/5_visualize.R \
        -s ${sim_dir}/stats.RDS \
        -o figures/simulation \
        -t

}

run_step5 () {

    echo ""
    echo "++ Running other scripts (tables, supplementary material) ++"
    echo ""

    Rscript code/other/get_demographics.R \
        -w ${proj_dir} \
        -o ${proj_dir}/tables

    Rscript code/other/get_results_table.R \
        -w ${proj_dir} \
        -o ${proj_dir}/tables

    Rscript code/other/get_excel.R \
        -w ${proj_dir} \
        -o ${proj_dir}/tables

    Rscript code/other/get_sim_results_table.R \
        -w ${proj_dir} \
        -o ${proj_dir}/tables
        
    Rscript code/other/get_supplementary.R

}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# get repo dir
proj_dir=$(git rev-parse --show-toplevel)

# get parameters for simulation
PARAM_IDX="1"  # index for simulation param file
param_path="${proj_dir}/code/3_Simulation/sim_param.txt"
n_param=$(cat ${param_path} | wc -l)

# get params from json
# https://unix.stackexchange.com/questions/459805/how-to-retrieve-values-from-json-object-using-awk-or-sed
n_simulation=$(cat ${proj_dir}/params.json | grep -o '"n_simulation":[^,}]*' | sed 's/.*: *\(.*\)/\1/')

# other paths
sim_dir="${proj_dir}/data/simulation"
feature_names_trunc="${proj_dir}/data/csv/feature_names_trunc.csv"

# ===============================
# ===== RUN PROJECT SCRIPTS =====
# ===============================

cd ${proj_dir}

for (( a=0; a<=5; a++ )); do

    stop_early "$a"

    if [ -z "$idx_to_run" ] || [ "$idx_to_run" -eq $a ]; then
        "run_step$a"
    fi

done
