#!/bin/bash
# 2 D with open boundary 1e7 runs up to sys 2^6. -1 means time() used for seed
#PBS -N BWS-2D-Open-N1e7
#PBS -m be
#PBS -q standard

#Adapt to your own folder structure if necessary
#cd ~/bws
#deterministic seed 42
./bws -S 42 -N 10000000 -L 6 -D 2 -B 0
