#!/bin/sh

RESULTSFILE=${1}
SRCDIR=${2}
BINDIR=${3}
SUITE=${4}

source ${SRCDIR}/Scripts/platforms.sh
touch parallel_test_output.txt
touch parallel_test_timings.txt


suiteFile=parallelUnitSuites
if [ ! -e ${suiteFile} ]; then
  suiteFile=${SRCDIR}/Scripts/parallelUnitSuites.txt
fi

# execute a particular suite if defined
if [ ! -z "${SUITE}" ]; then
    date
    printf "mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE}\n"
    date >> parallel_test_timings.txt
    printf "# mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE}\n" >> parallel_test_timings.txt
    mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE} >> parallel_test_output.txt
    date >> parallel_test_timings.txt
    printf "# ----------- \n"  >> parallel_test_timings.txt
    date
else
  for parallelUnitSuite in `cat ${suiteFile}`; do
    date
    printf "mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE}\n"
    date >> parallel_test_timings.txt
    printf "# mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE}\n" >> parallel_test_timings.txt
    mpirun -n 4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE} >> parallel_test_output.txt
    date >> parallel_test_timings.txt
    printf "# ----------- \n"  >> parallel_test_timings.txt
    date
  done
fi

#rm -f tmpresults_1.txt
#mpiexec -n 4 ${BINDIR}/plascom2_parallel_test -o tmpresults_1.txt
#
#@ i = 1
#while($i <= 8)
    #@ i += 1
    #if( -e tmpresults_1.txt ) then
        #@ i += 8;
    #else
        #sleep 30;
    #endif
#end
#
#if( -e tmpresults_1.txt ) then
  #cat tmpresults_1.txt >> ${RESULTSFILE}
#else
  #exit 1
#endif

exit 0


