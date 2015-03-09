#!/bin/bash
#$ -S /bin/bash
#$ -N deposition 
#$ -pe mpi 16
#$ -cwd
#$ -j y

module load openmpi-1.6.3-x86_64

python Deposition.py -i slowDepConfig.yml --debug --start 1002
