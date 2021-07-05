#!/bin/bash


BWS_DIR="/home/clustor/ma/g/gunsim/bws_code/BranchingWienerProcess"

for K in {1..10}
do
 qsub -v K=$K,SEED=$K run_histogram.sh
done 
