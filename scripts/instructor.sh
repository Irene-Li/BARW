#!/bin/bash

# 30/10/2018
# This is meant to be the replacement for job_dispatcher.sh and to be the main script for running new jobs of bww.
# Questions: b.walter16@imperial.ac.uk

# Change parameters below and then run with chmod +x instructor.sh && ./instructor
L=127
D=3
g=0 # Set g=1 for Sirpinski
ITER=256000000 # How many iterations
SEED=1005 # Starting seed
JOBS=256 # Subdivide in how many jobs?
C=100 # How many chunks per job

CODEDIR="/home/clustor/ma/g/gunsim/bws_code/BranchingWienerProcess"
SCRIPTDIR=${PWD}
OUTPUTDIR="/home/clustor/ma/g/gunsim/bws_output/"

SUBSCRIPTNAME="instructor_subscript.sh"

# ====================
# Below this line, it's not necessary to change things

# First, the number of iteratiosn is divided by the number of jobs

iterationperjob=$(echo " $ITER / $JOBS " | bc) 
version=$(${CODEDIR}/bws --version)

for i in $(seq 1 1 ${JOBS})
do
 #define joboutput
 datestring=$(date +%d-%m)
 JOBOUTPUT="${OUTPUTDIR}bws_run_${D}_${L}_${g}_${ITER}_${version}_${datestring}/bws_output_${SEED}.out"
 
 mkdir -p "$(dirname ${JOBOUTPUT})"
 echo "#!/bin/bash" > $SUBSCRIPTNAME
 echo "#PBS -N bws_${L}_${D}_${g}_$SEED" >> $SUBSCRIPTNAME
 echo "#PBS -m n" >> $SUBSCRIPTNAME
 echo "#PBS -q standard" >> $SUBSCRIPTNAME
 echo "cd ${CODEDIR}" >> $SUBSCRIPTNAME
 echo "./bws -N ${iterationperjob} -L $L -D $D -C $C --graph $g --seed ${SEED} > ${JOBOUTPUT}" >> $SUBSCRIPTNAME
 
 #qsub the job
 #chmod +x $SUBSCRIPTNAME 
 qsub $SUBSCRIPTNAME
 #change seed
 SEED=$((SEED + 10))
done


