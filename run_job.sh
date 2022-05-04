#!/bin/bash
# easier way to run bws.c
OUTPUTDIR="BARW/"
# Change parameters below and then run with chmod +x run_job.sh && ./run_job.sh
# sample: ./bws -L 257 -N 1000 -h 0.2 -p 0.25 -q 0.25 > BARW/data_h_0.4_p_0.25_q_0.25.txt 
SEED=1000 # Starting seed
L=256  # system size 2**n+1
N=1000 # number of realisations
hmin=0.3 # branching to hopping rate ratio
hmax=0.7 
hstep=0.1 

p=0.34 # prob. of hopping in the same direction as the previous hop 
q=0.33 # prob. of hopping to one of the two orthogonal directions.

# ====================
# Below this line, it's not necessary to change things
for i in $(seq ${hmin} ${hstep} ${hmax})
do
    #define joboutput
    JOBOUTPUT="${OUTPUTDIR}data_L${L}_N${N}_p${p}_q${q}_h${i}.dat"
    #mkdir -p "$(dirname ${JOBOUTPUT})"
    ./bws -L $L -N $N -h $i -p $p -q $q --seed ${SEED} > ${JOBOUTPUT}
    SEED=$((SEED + 10))
done