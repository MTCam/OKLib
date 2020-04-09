#!/bin/sh

outFile=${1}
srcDir=${2}
binDir=${3}

source ${srcDir}/Scripts/platforms.sh

ccrun -n 16 -t 5 ${binDir}/plascom2x -v 10 -p 

if [ ! -e ${outFile} ]; then
    printf "PlasCom2:RunsInParallel=" > ${outFile}
else
    printf "PlasCom2:RunsInParallel=" >> ${outFile}
fi

count=`grep Statistics PlasCom2Timing*txt | wc -l`
if [ ${count} == "1" ]; then
   printf "1\n" >> ${outFile}
else
   printf "0\n" >> ${outFile}
   exit 1
fi
exit 0
