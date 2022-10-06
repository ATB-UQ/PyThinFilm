#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -P m72 
#PBS -q gpuvolta
#PBS -N vac_dep
#PBS -l wd
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l walltime=48:00:00
#PBS -l mem=50GB
#PBS -l jobfs=1GB
export GMX_MAXBACKUP=-1
unset OMP_NUM_THREADS
module load openmpi/4.1.2
module load gromacs/2022-gpuvolta

${HOME}/.local/bin/pytf NM47_deposition.yml -n 1 
