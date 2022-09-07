#!/bin/bash
# easier way to run bws.c
OUTPUTDIR="BARW/sweep"
# Change parameters below and then run with chmod +x run_job.sh && ./run_job.sh
# sample: ./bws -L 257 -N 1000 -h 0.2 -p 0.25 -q 0.25 > BARW/data_h_0.4_p_0.25_q_0.25.txt 
SEED=1300 # Starting seed
L=257  # system size 2**n+1
N=10000 # number of realisations
hmin=0.2 # branching to hopping rate ratio
hmax=0.5
hstep=0.01

p=0.34 # prob. of hopping in the same direction as the previous hop 
q=0.33 # prob. of hopping to one of the two orthogonal directions.

a=0.75 # annihilation prob upon contact 

# ====================
# Below this line, it's not necessary to change things
for i in $(seq ${hmin} ${hstep} ${hmax})
do
    # define joboutput
    JOBOUTPUT="${OUTPUTDIR}data_L${L}_N${N}_p${p}_q${q}_h${i}_a${a}.dat"
    # mkdir -p "$(dirname ${JOBOUTPUT})"
    ./bws -L $L -N $N -h $i -p $p -q $q -a $a --seed ${SEED} > ${JOBOUTPUT}
    SEED=$((SEED + 10))
done