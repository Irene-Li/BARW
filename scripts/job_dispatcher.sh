#!/bin/bash
# This script provides a way to run various simulations for different values of L and h using the core-simulation file bws.c

# ENVIRONMENT PARAMETERS
BWS_DIR="/home/clustor/ma/g/gunsim/bws_code/BranchingWienerProcess"
OUTPUT_DIR="/home/clustor/ma/g/gunsim/bws_output/"

# SIMULATION PARAMETERS
# Write desired values of h
echo "0.1" > ${BWS_DIR}/hrange.list
# Dual logarithm of smallest and largest system size, i.e. LNMIN=5 marks start at 2^5-1 = 31
LNMIN=5
LNMAX=10
# Dimension
DIMENSION=5
# Boundary conditions (0 for open)
BOUNDARY=0
ITERATION=1000000
# Type of graph (0 for regular lattice, 1 for sierpinski, 2 for small-world network)
g=0


####################################
#for (( SIZE=2**${LNMIN} ; SIZE<=2**${LNMAX} ; SIZE*=2))
#do
#       for h in $(< ${BWS_DIR}/hrange.list)
#       do
#               qsub -v DIMENSION=${DIMENSION},SIZE=$((SIZE-1)),BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=${BOUNDARY},ITERATION=${ITERATION},h=$h,g=$g ${BWS_DIR}/scripts/bws_$#       done
#       let ITERATION=${ITERATION}*4
#done
##### SIRPINSKI
SIZEMIN=9
SIZEMAX=243
for (( SIZE=${SIZEMIN} ; SIZE<=${SIZEMAX} ; SIZE*=3))
do
        for h in $(< ${BWS_DIR}/hrange.list)
        do
                qsub -v DIMENSION=2,SIZE=${SIZE},BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=${BOUNDARY},ITERATION=${ITERATION},h=$h,g=1 ${BWS_DIR}/scripts/bws_wrapper.sh &
        done
        let ITERATION=${ITERATION}*4
done


rm ${BWS_DIR}/hrange.list
