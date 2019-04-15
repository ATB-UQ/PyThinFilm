#!/usr/bin/env bash

set -o errexit   #
set -o pipefail  # These lines cause bash to quit if any command has an error
set -o verbose
set -o nounset # error if variable is unset

name=$(basename $(pwd) )
datetime=$(date +"%F-%H-%M-%S")

mkdir -p log

nodes=3
cpus=$(bc<<<$nodes*28)

squeue --format "%Z %t" | grep $(pwd -P) | grep -v ' CG' || sbatch  <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${name}
#SBATCH --nodes=$nodes --ntasks-per-node=28
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=gpu
#SBATCH --qos=normal
#SBATCH --gres=gpu:tesla:2
#SBATCH --output log/${name}_${datetime}.stdout
#SBATCH --error log/${name}_${datetime}.stderr
set -o errexit
set -o pipefail
set -o verbose
#load modules here
module use /scratch/scmb/uqtlee10/modulefiles
module load gromacs/2018.3-TLmod
python deposition.py --name $name -i parameters.yml --debug --max-cores $cpus >> log/${name}_dep.log 2>&1
EOF
#84 56 28
