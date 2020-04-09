#!/bin/bash

OutFile=${1}
SRCDIR=${2}
BINDIR=${3}
TmpOut=${OutFile}_tmp.txt

HOSTNAME=`hostname -s`
RUNCOMMAND=""
BATCHCOMMAND=""

set echo 
rm -rf plascom2x_serial.out

if [[ ${HOSTNAME} = "vulcan"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "cab"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "quartz"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "syrah"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "golub"* ]]; then
  RUNCOMMAND="ccrun -n 1 -t 2 --"
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "taub"* ]]; then
  RUNCOMMAND="ccrun -n 1 -t 2 --"
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "knl"* ]]; then
  RUNCOMMAND="srun -n 1 -p development -t 2 -A TG-CTS090004 ibrun "
  ${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
elif [[ ${HOSTNAME} = "titan"* ]]; then
  RUNCOMMAND="aprun -n1"
  BATCHCOMMAND="qsub -I -A csc188 -l nodes=1:ppn=1,walltime=00:05:00 -x"
  set echo 
  rm -rf plascom2x_serial.out
  rm -f plascom2x_serial_batch.sh
  printf "#!/bin/bash\n" > plascom2x_serial_batch.sh
  printf "cd \${PBS_O_WORKDIR}n" >> plascom2x_serial_batch.sh
  printf "${RUNCOMMAND} ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out\n" >> plascom2x_serial_batch.sh
  chmod +x plascom2x_serial_batch.sh
  ${BATCHCOMMAND} ./plascom2x_serial_batch.sh
else
  ${BINDIR}/plascom2x -v 3 -p >& plascom2x_serial.out
fi

# Test the serial example program output
printf "PlasCom2:RunsInSerial=" > ${TmpOut}
if [ ! -e "PlasCom2Timing_000001.txt" ]; then
   echo "Couldn't find plascom2x profiling result."
   printf "0\n" >> ${TmpOut}
else
  printf "1\n" >> ${TmpOut}
#  printf "PlasCom2:Runs=1\n" >> ${TmpOut}
fi
STEST=`cat plascom2x_serial.out | grep Hello`
printf "PlasCom2:CanCallFortran=" >> ${TmpOut}
if [ -z "${STEST}" ]; then
  printf "0\n" >> ${TmpOut}
else
  printf "1\n" >> ${TmpOut}
fi

rm -rf PlasCom2Timing_000001.txt
if [ ! -e ${OutFile} ]; then
    cat ${TmpOut} > ${OutFile}
else
    cat ${TmpOut} >> ${OutFile}
fi
rm -f ${TmpOut}
rm -f plascom2x_serial.out
exit 0
