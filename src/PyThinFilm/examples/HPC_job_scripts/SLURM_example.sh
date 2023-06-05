#!/bin/bash --login
#SBATCH --job-name=PyThinFilm_Example
#SBATCH --account=m72
#SBATCH --partition=work
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:10:00

module load gromacs/2021.4
module load hpc-python-collection/2022.11-py3.9.15
export OMP_NUM_THREADS=8

# Temporal workaround for avoiding Slingshot issues on shared nodes:
export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)

export GMX_MAXBACKUP=-1
export MPICH_OFI_STARTUP_CONNECT=1
export NRANK=$(($SLURM_NTASKS/$OMP_NUM_THREADS))

pytf SLURM_solvent_evaporation.yml
