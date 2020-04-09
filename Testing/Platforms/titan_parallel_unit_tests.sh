#!/bin/sh

RESULTSFILE=${1}
SRCDIR=${2}
BINDIR=${3}
SUITE=${4}

source ${SRCDIR}/Scripts/titan.sh

if [ ! -z "${SUITE}" ]; then
    rm -f ${SUITE}_batch.sh
    printf "#!/bin/sh\n" > ${SUITE}_batch.sh
    printf "cd \${PBS_O_WORKDIR}\n" >> ${SUITE}_batch.sh
    printf "${RUNCOMMAND} -n4 ${BINDIR}/plascom2_parallel_test -n ${SUITE} -o ${RESULTSFILE}\n" >> ${SUITE}_batch.sh
    chmod +x ./${SUITE}_batch.sh
    ${BATCHCOMMAND} ./${SUITE}_batch.sh
else
  for parallelUnitSuite in `cat ${SRCDIR}/Scripts/parallelunitsuites.txt`; do
    rm -f ${parallelUnitSuite}_batch.sh
    printf "#!/bin/sh\n" > ${parallelUnitSuite}_batch.sh
    printf "cd \${PBS_O_WORKDIR}\n" >> ${parallelUnitSuite}_batch.sh
    printf "${RUNCOMMAND} -n4 ${BINDIR}/plascom2_parallel_test -n ${parallelUnitSuite} -o ${RESULTSFILE}\n" >> ${parallelUnitSuite}_batch.sh
    chmod +x ./${parallelUnitSuite}_batch.sh
    ${BATCHCOMMAND} ./${parallelUnitSuite}_batch.sh
  done
fi
exit 0
