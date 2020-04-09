#!/bin/bash

# This is a testing stub - it is not really a test. There are no 
# real general platform-specific tests right now.

if [ ! -e ${1} ]; then
    printf "TestStubWorks=1\n" > ${1}
else 
    printf "TestStubWorks=1\n" >> ${1}
fi

exit 0
