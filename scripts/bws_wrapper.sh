#!/bin/bash
#PBS -N bws_run_${DIMENSION}_${BOUNDARY}_${SIZE}_$h
#PBS -q standard
#PBS -m n
# For variable name explanation see job_dispatcher.sh

echo `date` >> ${OUTPUT_DIR}/log.dat
echo "./bws -D ${DIMENSION} --BCs ${BOUNDARY} -L ${SIZE} -N ${ITERATION} -h $h --graph $g > ${OUTPUT_DIR}/bws_$PBS_JOBID.out" >> ${OUTPUT_DIR}/log.dat
cd ${BWS_DIR}

VERSION="$(./bws --version)"
echo "${DIMENSION} ${BOUNDARY} ${VERSION}" >> ${OUTPUT_DIR}/log.dat
./bws -D ${DIMENSION} --BCs ${BOUNDARY} -L ${SIZE} -N ${ITERATION} -h $h --graph $g > ${OUTPUT_DIR}/bws_${VERSION}_$PBS_JOBID.out 

echo "---------------------" >> ${OUTPUT_DIR}/log.dat

