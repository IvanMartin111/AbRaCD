#!/bin/bash

###############################################################
# CONFIGURATION
###############################################################

NP=10
REPS=3

NLOOPS=20000
LOCOBEST=200

ENERGY="korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin"

###############################################################
# FUNCTION
###############################################################

run_cdr () {

    local CDR=$1
    local FLAG=$2
    local KINK=$3
    local DVAL=$4

    local INPUT="AB_loop_RCD_${CDR}_XRayBench.txt"

    echo
    echo "=========================================="
    echo "Running $CDR"
    echo "=========================================="
    echo

    ###########################################################
    # RCD
    ###########################################################

    for REP in $(seq 1 $REPS)
    do

        local OUT="20K_${CDR}_NBench_map_rep${REP}"

        echo "Launching $OUT"

        mpirun -np $NP rcd_mpi \
            $INPUT \
            --ABmap CDRS.bin \
            $FLAG \
            --nokink_merge $KINK \
            -r \
            -e 56 \
            --energy_file $ENERGY \
            --bench \
            -n $NLOOPS \
            --loco_best $LOCOBEST \
            -x dunbrack.bin \
            -d $DVAL \
            --linear \
            -t 0.99 \
            -o $OUT

    done

    ###########################################################s
    # REFINEMENT
    ###########################################################

    for REP in $(seq 1 $REPS)
    do

        DIR="20K_${CDR}_NBench_map_rep${REP}"

        echo "Refining $DIR"

        cd "$DIR" || exit

        for IVF in $(seq 1 $REPS)
        do

            refine_loops.pl \
                ../$INPUT \
                native \
                _native_20K_${CDR}_NBench_map_ivf_rep${IVF} \
                15 \
                "--energy ref2015"

        done

        cd ..

    done

}

###############################################################
# MAIN
###############################################################

run_cdr H1 --H1 0.4 0.2

run_cdr H2 --H2 0.9 0.5

run_cdr L1 --L1 0.9 0.5

run_cdr L2 --L2 0.9 0.5

run_cdr L3 --L3 0.9 0.5

echo
echo "=========================================="
echo "All calculations finished."
echo "=========================================="
echo
