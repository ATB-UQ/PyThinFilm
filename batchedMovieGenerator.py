from MovieGenerator import getSortedRunDirList
from Deposition import PROJECT_DIR
from os.path import join, exists, basename
import os
import yaml
import argparse
import subprocess
from jinja2 import Template

JOBS_DIRECTORY = "qsub"
MAX_N_CORES = 100

JOB_TEMPLATE = '''#!/bin/bash
{%- if n_cores %}
#$ -pe mpi {{n_cores}}
{%- endif %}
#$ -S /bin/bash
#$ -N Mov{{batchStr_escaped}}
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

python MovieGenerator.py -i {{config_file}} -b {{batchStr}}
'''

def runBatchMovie(runConfig, args):
    
    workdir = join(PROJECT_DIR, runConfig["work_directory"])
    
    jobDir = join(workdir, JOBS_DIRECTORY)
    if not exists(jobDir):
        os.mkdir(jobDir)
    fullRunList = getSortedRunDirList(workdir, args.batch)
    numberOfJobs = len(fullRunList)

    
    # integer division on purpose
    n_cores = runConfig['movies']['n_cores']
    jobsPerCore = (numberOfJobs * n_cores)/(MAX_N_CORES + 1) + 1

    for dirname in chunks(fullRunList, jobsPerCore):
        i = int(basename(dirname[0]))
        print i
        start, end = i, i + jobsPerCore - 1
        t = Template(JOB_TEMPLATE)
        jobStr = t.render( config_file=args.input, batchStr="{0}:{1}".format(start, end), batchStr_escaped="{0}-{1}".format(start, end), n_cores=n_cores )
        
        jobPath = join(jobDir, "BM_{0}-{1}.sh".format(start, end))
        with open(jobPath, "w") as fh:
            fh.write(jobStr)
        subprocess.Popen("qsub {0}".format(jobPath).split()).wait()


# Source : http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def parseCommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-b', '--batch', dest='batch')
    
    args = parser.parse_args()
    runConfig = yaml.load(open(args.input))
    
    runBatchMovie(runConfig, args)
    
parseCommandline()
