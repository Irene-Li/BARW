#!/bin/bash

BWS_DIR="/home/clustor/ma/g/gunsim/bws_code/BranchingWienerProcess"
OUTPUT_DIR="/home/clustor/ma/g/gunsim/bws_output/"


#qsub -v DIMENSION=2,SIZE=27,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=1000000,h=0.1,g=1 ${BWS_DIR}/scripts/bws_wrapper.sh &
#qsub -v DIMENSION=2,SIZE=81,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=1000000,h=0.1,g=1 ${BWS_DIR}/scripts/bws_wrapper.sh &
#qsub -v DIMENSION=2,SIZE=243,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=1000000,h=0.1,g=1 ${BWS_DIR}/scripts/bws_wrapper.sh &
#qsub -v DIMENSION=2,SIZE=729,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=4000000,h=0.1,g=1 ${BWS_DIR}/scripts/bws_wrapper.sh &
#qsub -v DIMENSION=2,SIZE=255,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=64000000,h=0.1,g=0,SEED=10 ${BWS_DIR}/scripts/bws_wrapper.sh &
qsub -v DIMENSION=5,SIZE=63,BWS_DIR="${BWS_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BOUNDARY=0,ITERATION=64000000,h=0.1,g=0,SEED=10 ${BWS_DIR}/scripts/bws_wrapper.sh &

