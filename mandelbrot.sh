#!/bin/bash

#PBS -l select=2:ncpus=2:mem=2gb

#PBS -l walltime=0:01:00

#PBS -q short_cpuQ

module load mpich-3.2
cd ./HPC4DS-Project
mpirun.actual -n 2 ./mandelbrot 1 100 #number of zooms, number of threads