#!/bin/bash

#
# this is a script to run a program and return 0 if the exit code
# is equal to the first argument.
#

# save the expected exit code
expected_code=$1
shift

# run the command
$*

# grab  the exit code
exit_code=$?

if [ $exit_code -eq $expected_code ]; then
    exit 0
else
    echo "exit code for $1 was $exit_code, expected $expected_code"
    exit 1
fi
