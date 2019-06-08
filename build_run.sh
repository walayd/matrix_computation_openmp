#!/usr/bin/env bash

# build solution
mpicc -std=c99 -Wall -lm -fopenmp -o out  larabi.c

# run with some dummy data

# mpirun --oversubscribe -np 2 ./out  ./IFI-RIF-Projet2/data/a_4 ./IFI-RIF-Projet2/data/b_4
# mpirun --oversubscribe -np 4 ./out  ./IFI-RIF-Projet2/data/a_4 ./IFI-RIF-Projet2/data/b_4

# mpirun --oversubscribe -np 2 ./out  ./IFI-RIF-Projet2/data/a_8 ./IFI-RIF-Projet2/data/b_8
mpirun --oversubscribe -np 4 ./out  ./IFI-RIF-Projet2/data/a_8 ./IFI-RIF-Projet2/data/b_8
# mpirun --oversubscribe -np 8 ./out  ./IFI-RIF-Projet2/data/a_8 ./IFI-RIF-Projet2/data/b_8

# cd ./IFI-RIF-Projet2
# python3 ./evaluate.py
# cd ..

# clear solution
rm ./out