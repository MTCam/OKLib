#!/bin/sh

outFile=${1}
srcDir=${2} 
binDir=${3}

source ${srcDir}/Scripts/platforms.sh

suiteFile=serialUnitSuites
if [ ! -e ${suiteFile} ]; then
  suiteFile=${srcDir}/Scripts/serialUnitSuites.txt
fi

for serialUnitSuite in `cat ${suiteFile}`; do
  if [ -z "${LSB_HOSTS}" ]; then
    printf "lalloc 1 -q pbatch lrun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
    lalloc 1 -q pbatch lrun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}
  else
    printf "lrun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
    lrun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}
  fi
done

printf "Done running serial Suites\n"
#srun -n 1  -ppdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config 
#srun -n 2  -ppdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 4  -ppdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 8  -ppdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 16 -ppdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config


