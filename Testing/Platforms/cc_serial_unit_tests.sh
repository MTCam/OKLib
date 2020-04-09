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
  printf "ccrun -n 1 ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
  ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}
done

printf "Done running serial Suites\n"


