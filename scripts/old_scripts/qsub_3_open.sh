#!/bin/bash
#2 D with open boundary 100,000 runs up to sys 2^6. -1 means time() used for seed
#PBS -N BWSStandard5
#PBS -m be
#PBS -q standard
cd ~/bws
./bws -1 100000 9 3 0
