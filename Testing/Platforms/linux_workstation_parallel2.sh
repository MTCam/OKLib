#!/bin/bash

RESULTSFILE=${1}
SRCDIR=${2}
BINDIR=${3}

source ${SRCDIR}/Scripts/platforms.sh

rm -f PlasCom2Timing_000004.txt

mpiexec -n 4 ${BINDIR}/plascom2x -v 10 -p 

if [ ! -e ${RESULTSFILE} ]; then
    printf "PlasCom2:RunsInParallel=" > ${RESULTSFILE}
else
    printf "PlasCom2:RunsInParallel=" >> ${RESULTSFILE}
fi

err=0

if [ ! -e PlasCom2Timing_000004.txt ]; then
    printf "PlasCom2:RunsInParallel:Error: Failed to produce expected output timing file.\n"
    err=1
fi

count=`grep Statistics PlasCom2Timing_000004.txt | wc -l`

if [ ${count} == "1" ]; then
   printf "1\n" >> ${RESULTSFILE}
else
   printf "PlasCom2:RunsInParallel:Error: Unexpected count = ${count}\n"
   printf "0\n" >> ${RESULTSFILE}
   err=1
fi

exit ${err}
