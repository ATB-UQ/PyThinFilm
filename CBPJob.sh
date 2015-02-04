#!/bin/bash
#$ -S /bin/bash
#$ -N testmpi 
#$ -pe mpi 16
#$ -cwd
#$ -j y

module load openmpi-1.6.3-x86_64

python CBPDeposition.py
