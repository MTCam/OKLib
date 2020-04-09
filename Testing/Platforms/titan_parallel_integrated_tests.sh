#!/bin/sh

outFile=${1}
srcDir=${2}
binDir=${3}

source ${srcDir}/Scripts/titan.sh

rm -f first_integrated_test.sh
printf "#!/bin/sh\n" > first_integrated_test.sh
printf "cd \${PBS_O_WORKDIR}\n" >> first_integrated_test.sh
printf "${RUNCOMMAND} -n16 ${binDir}/plascom2x -v 10 -p\n" >> first_integrated_test.sh
chmod +x ./first_integrated_test.sh
${BATCHCOMMAND} ./first_integrated_test.sh


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
