#!/bin/bash
#PBS -N bws_127_3_0_3555
#PBS -m n
#PBS -q standard
cd /home/clustor/ma/g/gunsim/bws_code/BranchingWienerProcess
./bws -N 1000000 -L 127 -D 3 -C 100 --graph 0 --seed 3555 > /home/clustor/ma/g/gunsim/bws_output/bws_run_3_127_0_256000000_v1-120-g9366_31-10/bws_output_3555.out
