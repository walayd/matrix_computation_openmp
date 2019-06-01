#!/usr/bin/env bash

# build solution
mpicc -Wall -lm -fopenmp -o out  larabi.c

# run with some dummy data
mpirun  ./out  ./dummy/data ./dummy/data2 --oversubscribe

# clear solution
rm ./out