#!/usr/bin/env bash

set -o errexit   #
set -o pipefail  # These lines cause bash to quit if any command has an error
set -o verbose
set -o nounset # error if variable is unset

name=$(basename $(pwd) )
datetime=$(date +"%F-%H-%M-%S")

mkdir -p log

ncpus=24
ngpus=8

qstat -f |tr -d " \t\n\r" | grep $(pwd) > /dev/null || qsub  <<EOF
#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -P m72
#PBS -q gpu
#PBS -N ${name}
#PBS -l wd
#PBS -l ncpus=${ncpus}
#PBS -l ngpus=${ngpus}
#PBS -l walltime=47:30:00
#PBS -l mem=24GB
#PBS -l jobfs=100MB
#PBS -o log/${name}_${datetime}.stdout
#PBS -e log/${name}_${datetime}.stderr
set -o errexit
set -o pipefail
set -o verbose
#load modules here
module load gromacs/2018.3-gpu
python deposition.py --name $name -i parameters.yml --debug --max-cores ${ncpus} >> log/${name}_dep.log 2>&1
EOF
