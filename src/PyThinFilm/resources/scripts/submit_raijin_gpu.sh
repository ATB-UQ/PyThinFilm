#!/usr/bin/env bash

set -o errexit   #
set -o pipefail  # These lines cause bash to quit if any command has an error
set -o verbose

name=$(basename $(pwd) )
datetime=$(date +"%F-%H-%M-%S")

mkdir -p log

ncpus=6
ngpus=2

if [ "$1" == "force" ]
then
    qstat="echo x"
else
    qstat="qstat -f"
fi

$qstat | tr -d ' \t\n\r' | grep $(pwd) > /dev/null || qsub  <<EOF
#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -P m72
#PBS -q gpu
#PBS -N ${name}
#PBS -l wd
#PBS -l ncpus=${ncpus}
#PBS -l ngpus=${ngpus}
#PBS -l walltime=12:00:00
#PBS -l mem=24GB
#PBS -l jobfs=100MB
#PBS -o log/${name}_${datetime}.stdout
#PBS -e log/${name}_${datetime}.stderr
set -o errexit
set -o pipefail
set -o verbose

term_handler()
{
    echo "Restarting terminated job..."
    ./submit_raijin_gpu.sh force
    echo "Restarted terminated job."
    exit -1
}

# uncomment to enable term handler
# trap 'term_handler' SIGTERM

export GMX_MAXBACKUP=-1

#load modules here
module load gromacs/2018.3-gpu
python deposition.py --name $name -i parameters.yml --debug --max-cores ${ncpus}
EOF
