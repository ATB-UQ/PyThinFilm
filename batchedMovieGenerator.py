from MovieGenerator import getSortedRunDirList
from Deposition import PROJECT_DIR
from os.path import join, exists
import os
import yaml
import argparse
import subprocess

JOBS_DIRECTORY = "qsub"
MAX_N_CORES = 100

JOB_TEMPLATE = '''#!/bin/bash
#$ -S /bin/bash
#$ -N batchMovie 
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

python MovieGenerator.py -i {0} -b {1}
'''

def runBatchMovie(runConfig, args):
    
    workdir = join(PROJECT_DIR, runConfig["work_directory"])
    
    jobDir = join(workdir, JOBS_DIRECTORY)
    if not exists(jobDir):
        os.mkdir(jobDir)
    fullRunList = getSortedRunDirList(workdir, "")
    numberOfJobs = len(fullRunList)
    
    # integer division on purpose
    jobsPerCore = numberOfJobs/(MAX_N_CORES + 1) + 1
    
    i = 1
    while i < len(fullRunList):
        start, end = i, i + jobsPerCore - 1
        jobStr = JOB_TEMPLATE.format(args.input, "{0}:{1}".format(start, end))
        
        jobPath = join(jobDir, "BM_{0}-{1}.sh".format(start, end))
        with open(jobPath, "w") as fh:
            fh.write(jobStr)
        subprocess.Popen("qsub {0}".format(jobPath).split()).wait()
        
        i += jobsPerCore
    

def parseCommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    
    args = parser.parse_args()
    runConfig = yaml.load(open(args.input))
    
    runBatchMovie(runConfig, args)
    
parseCommandline()
