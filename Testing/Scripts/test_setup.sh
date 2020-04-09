#!/bin/sh

OUTFILE=${1}
SRCDIR=${2}
BINDIR=${3}

date > setupoutput
pwd >> setupoutput
ls -lat * >> setupoutput
rm -rf *.out
rm -rf *.txt
rm -rf *.h5
rm -rf *.xdmf
rm -rf *.txt
ls -lat * >> setupoutput
printf "Removed all previous testing results.\n" >> setupoutput
cp -r ${SRCDIR}/Data/* .
printf "Staged all testing data to testing area.\n" >> setupoutput
ls -lat * >> setupoutput
mv setupoutput ${OUTFILE}
exit 0

