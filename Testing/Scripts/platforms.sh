#!/bin/sh

HOSTNAME=`hostname -s`
PRUNCOMMAND=""
RUNCOMMAND=""
PRUNARGS=""

if [[ "${HOSTNAME}" == "vulcan"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  PRUNCOMMAND="srun"
  PRUNARGS="-ppdebug"
fi
if [[ "${HOSTNAME}" == "cab"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  PRUNCOMMAND="srun"
  PRUNARGS="-ppdebug"
fi
if [[ "${HOSTNAME}" == "quartz"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  PRUNCOMMAND="srun"
  PRUNARGS="-ppdebug"
fi
if [[ "${HOSTNAME}" == "syrah"* ]]; then
  RUNCOMMAND="srun -n 1 -ppdebug -t 2 "
  PRUNCOMMAND="srun"
  PRUNARGS="-ppdebug"
fi
if [[ "${HOSTNAME}" == "golub"* ]]; then
  RUNCOMMAND="ccrun -n 1 -t 2 --"
  PRUNCOMMAND="ccrun"
  PRUNARGS="--"
fi
if [[ "${HOSTNAME}" == "taub"* ]]; then
  RUNCOMMAND="ccrun -n 1 -t 2 --"
  PRUNCOMMAND="ccrun"
  PRUNARGS="--"
fi
if [[ "${HOSTNAME}" == "knl"* ]]; then
  RUNCOMMAND="srun -n 1 -p development -t 2 -A TG-CTS090004 ibrun "
  PRUNCOMMAND="srun"
  PRUNARGS="-p development -A TG-CTS090004 ibrun"
fi
export PRUNCOMMAND PRUNARGS

