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

# execute a particular suite if defined
if [ ! -z "${SUITE}" ]; then
    date
    printf "bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt\n"
    date >> parallel_test_timings.txt
    printf "#bsub -n 20  -q pdebug -G guests -Is -XF -W60mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
    time bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt
    date >> parallel_test_timings.txt
    printf "# --------- \n" >> parallel_test_timings.txt
    date
else
  for parallelUnitSuite in `cat ${suiteFile}`; do
    date
    printf "bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parllelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n"
    date >> parallel_test_timings.txt
    printf "# bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parllelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
    time bsub -n 20 -q pdebug -G guests -Is -XF -W60 mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parllelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt
    date >> parallel_test_timings.txt
    printf "# --------- \n" >> parallel_test_timings.txt
    date
  done
fi
exit 0
