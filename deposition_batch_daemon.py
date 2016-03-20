import yaml
import argparse
import os
from os.path import dirname, exists, abspath, basename
import subprocess
import logging

CLUSTER_TEMPLATE_NAME = "{name}_run_script.sh.template"
BATCH_RUNS_DIR = "batch_runs"
PBS_GET_ALL_MY_JOB_NAMES = 'if [ -n "$(qselect -u $USER)" ] ; then qselect -u $USER | xargs qstat -f | grep Job_Name ; fi'

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')

if not exists(BATCH_RUNS_DIR):
    os.makedirs(BATCH_RUNS_DIR)

def run(args, cwd_opt=None):
    return subprocess.Popen(args, cwd=cwd_opt).wait()

def run_with_shell_stdout_piped(args):
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
    stdout, _ = proc.communicate()
    return stdout

def setup_batch_config(config_file):
    config = yaml.load(open(config_file))
    return batch_config_file

def populate_run_script_template(template_name, config, batch_config_file):
    cluster_template = open(CLUSTER_TEMPLATE_NAME.format(name=template_name)).read()
    return run_script

def submit_to_queue(run_script):
    run(["qsub", run_script], cwd_opt=dirname(abspath(run_script)))

def job_in_queue(run_script):
    qstat = run_with_shell_stdout_piped(PBS_GET_ALL_MY_JOB_NAMES)
    return basename(run_script) in qstat

def get_running_job(config_file):
    return "test.sh"
    config = yaml.load(open(config_file))
    # TODO: check if batched job exists and if so whether last one is running
    job_in_queue(run_script)

def run_batched_deposition(master_config_file, cluster_template_name):
    running_job = get_running_job(master_config_file, )
    if not running_job:
        master_config, batch_config_file = setup_batch_config(args.input)
        run_script = populate_run_script_template(cluster_template_name, master_config, batch_config_file)
        submit_to_queue(run_script)
    else:
        logging.info("A batch job is already running: {0}".format(running_job))

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--master_config')
    parser.add_argument('-c', '--cluster_template')
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args = parse_command_line()
    run_batched_deposition(args.master_config, args.cluster_template)