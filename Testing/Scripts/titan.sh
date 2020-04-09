#!/bin/sh

ACCOUNT="csc188"
BATCHCOMMAND="qsub -I -A ${ACCOUNT} -q debug -l nodes=1,walltime=01:00:00 -x"
RUNCOMMAND="aprun"
