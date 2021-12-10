#!/bin/bash

#PBS -l select=16:ncpus=8:mem=2gb

#PBS -l walltime=0:01:00

#PBS -q short_cpuQ

module load mpich-3.2
cd ./HPC4DS-Project

mpirun.actual -n 128 ./mandelbrot -z 8 -t 100 #number of zooms, number of threads

