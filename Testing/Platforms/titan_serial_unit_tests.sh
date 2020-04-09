#!/bin/sh

outFile=${1}
srcDir=${2} 
binDir=${3}

source ${srcDir}/Scripts/titan.sh

suiteFile=serialUnitSuites
if [ ! -e ${suiteFile} ]; then
  suiteFile=${srcDir}/Scripts/serialUnitSuites.txt
fi

for serialUnitSuite in `cat ${suiteFile}`; do
  #printf "srun -n 1 -ppdebug -t 5 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
  #srun -n 1 -ppdebug -t 5 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}

  rm -f ${serialUnitSuite}_batch.sh
  printf "#!/bin/sh\n" > ${serialUnitSuite}_batch.sh
  printf "cd \${PBS_O_WORKDIR}\n" >> ${serialUnitSuite}_batch.sh
  printf "${RUNCOMMAND} -n1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n" >> ${serialUnitSuite}_batch.sh
  ${BATCHCOMMAND} ./${serialUnitSuite}_batch.sh
done

