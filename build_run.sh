#!/usr/bin/env bash

# build solution
mpicc -Wall -lm -fopenmp -o out  larabi.c

# run with some dummy data
# mpirun --oversubscribe -np 1 ./out  ./dummy/data ./dummy/data2
mpirun --oversubscribe ./out  ./dummy/data ./dummy/data2

# clear solution
rm ./out