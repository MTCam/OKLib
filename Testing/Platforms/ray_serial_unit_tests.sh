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
  printf "bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
  bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}
done

printf "Done running serial Suites\n"
#srun -n 1  -npdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config 
#srun -n 2  -npdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 4  -npdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 8  -npdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config
#srun -n 16 -npdebug -t 5 ${binDir}/advect1d -p 2 -c ./advect1d.config


