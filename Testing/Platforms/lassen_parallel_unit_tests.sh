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
    printf "lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt\n"
    date >> parallel_test_timings.txt
    if [ -z "${LSB_HOSTS}" ]; then
      printf "# lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
      time lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt
    else
      printf "# lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
      time lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt
    fi
    date >> parallel_test_timings.txt
    printf "# --------- \n" >> parallel_test_timings.txt
    date
else
  for parallelUnitSuite in `cat ${suiteFile}`; do
    date
    if [ -z "${LSB_HOSTS}" ]; then
      printf "lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n"
      date >> parallel_test_timings.txt
      printf "# lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
      time lalloc 1 lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt
    else
      printf "lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n"
      date >> parallel_test_timings.txt
      printf "# lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt\n" >> parallel_test_timings.txt
      time lrun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt
    fi
    date >> parallel_test_timings.txt
    printf "# --------- \n" >> parallel_test_timings.txt
    date
  done
fi
exit 0