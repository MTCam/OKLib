#!/bin/sh

outFile=${1}
srcDir=${2} 
binDir=${3}

source ${srcDir}/Scripts/platforms.sh

suiteFile=serialUnitSuites
if [ ! -e ${suiteFile} ]; then
  suiteFile=${SRCDIR}/Scripts/serialUnitSuites.txt
fi

for serialUnitSuite in `cat ${suiteFile}`; do
  printf "${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}\n"
  ${binDir}/plascom2_test -n ${serialUnitSuite} -o ${outFile}
done



