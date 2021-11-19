#!/bin/bash

if [ "$#" -eq 1 ]; then
    module load mpich-3.2
    fbname=$(basename "$1" .c)
    fpath=$(dirname "$1")

    set -x #To echo
    mpicc -std=c99 -g -Wall -o $fpath/$fbname $1 -lm
else
    echo "Utilizzo: ./compile.sh <path_to_c_file>"
fi