#!/bin/bash
#PBS -P m72
#PBS -q normal
#PBS -l ncpus=32
#PBS -l walltime=48:00:00
#PBS -l mem=16Gb
#PBS -l wd


module load python
module load gromacs/4.0.7 
module load openmpi/1.6.3

python Deposition.py -i CBP-10nm-CBP-10nm-450K-slow.yml --start 2190
