#!/bin/csh

setenv pwd ${PWD}
make clean
make
root.exe -q -l -b ${pwd}/run.C 
