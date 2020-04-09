#!/bin/sh

RESULTSFILE=${1}
SRCDIR=${2}
BINDIR=${3}
SUITE=${4}

source ${SRCDIR}/Scripts/platforms.sh
touch parallel_test_timings.txt
touch parallel_test_output.txt

suiteFile=parallelUnitSuites
if [ ! -e ${suiteFile} ]; then
  suiteFile=${SRCDIR}/Scripts/parallelUnitSuites.txt
fi

if [ ! -z "${SUITE}" ]; then
  date
  printf "ccrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE}\n"
  date >> parallel_test_timings.txt
  printf "# ccrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE}\n" >> parallel_test_timings.txt
  time ccrun -n 4 -- ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt 
  date >> parallel_test_timings.txt
  printf "# --------- \n" >> parallel_test_timings.txt
  date
else
  for parallelUnitSuite in `cat ${suiteFile}`; do
    date
    printf "ccrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE}\n"
    date >> parallel_test_timings.txt
    printf "# ccrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE}\n" >> parallel_test_timings.txt
    time ccrun -n 4 -- ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt 
    date >> parallel_test_timings.txt
    printf "# --------- \n" >> parallel_test_timings.txt
    date
  done
fi

exit 0
