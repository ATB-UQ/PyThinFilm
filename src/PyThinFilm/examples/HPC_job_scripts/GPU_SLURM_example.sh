#!/bin/bash --login
#SBATCH --job-name=PyThinFilm_Example
#SBATCH --partition=gpu
#SBATCH --nodes=1              #1 nodes in this example
#SBATCH --ntasks-per-node=1    #1 tasks for the 1 GPUs in this job
#SBATCH --gpus-per-node=1      #1 GPUs in this job
#SBATCH --sockets-per-node=1   #Use the 1 slurm-sockets in this job
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --account=m72-gpu

module load PrgEnv-cray
module load rocm craype-accel-amd-gfx90a

module load gcc/12.1.0
module load gromacs-amd-gfx90a/2022.3.amd1_174
module load hpc-python-collection/2022.11-py3.9.15

unset OMP_NUM_THREADS
export GMX_MAXBACKUP=-1

pytf GPU_solvent_evaporation.yml

# mpi_template : "srun -l -u -c 8"
# mdrun_template: "{GMX_EXEC} mdrun -s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final} -nb gpu -bonded gpu -pin on -update gpu -ntomp 8 -ntmpi 1"
# --> srun -l -u -c 8 gmx mdrun -s ... -nb gpu -bonded gpu -pin on -update gpu -ntomp 8 -ntmpi 1
