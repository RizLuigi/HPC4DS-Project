#!/bin/bash

#PBS -l select=2:ncpus=2:mem=2gb

#PBS -l walltime=0:02:00

#PBS -q short_cpuQ

module load mpich-3.2
cd ./mandelbrot
mpirun.actual -n 4 ./mandelbrot 16 100 #number of zooms, number of threads