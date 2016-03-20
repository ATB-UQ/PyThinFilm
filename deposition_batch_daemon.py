import yaml
import argparse
import os
from os.path import dirname, exists, abspath, basename, join, isdir 
import subprocess
import logging
import glob

CLUSTER_TEMPLATE_NAME = "{name}_run_script.sh.template"
BATCH_RUN_SCRIPT_NAME = "deposition_batch_{start}_{end}.sh"

BATCH_RUNS_DIR = "batch_runs"
BATCH_RUN_DIR_NAME = "from_{start}_to_{end}"
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
    return
    qstat = run_with_shell_stdout_piped(PBS_GET_ALL_MY_JOB_NAMES)
    if basename(run_script) in qstat:
        return basename(run_script)

def is_being_processed(config):
    run_dirs = [join(config["work_directory"], f) for f in os.listdir(config["work_directory"]) if isdir(join(config["work_directory"], f))]
    last_run_dir = sorted(run_dirs)[-1] if run_dirs else None
    logging.info("Last run dir: {0}".format(last_run_dir))
    last_deposition_steps = sorted(config["deposition_steps"], key=lambda x:x["first_sim_id"])[-1]
    logging.info("Last deposition step in master config: {0}".format(last_deposition_steps["description"]))
    # check whether last_run_dir is within the scope of the last deposition step
    last_deposition_being_processed = last_deposition_steps["first_sim_id"] <= int(basename(last_run_dir)) < last_deposition_steps["last_sim_id"]
    if not last_deposition_being_processed:
        if int(basename(last_run_dir)) + 1 == last_deposition_steps["first_sim_id"]:
            # can initialise last deposition step
            return False
        else:
            # something doesn't match up
            raise Exception("Cannot proceed as run dirs don't match last deposition step in master config")
    return True

def get_last_batch_run_dir():
    batch_run_dirs = sorted([join(BATCH_RUNS_DIR, f) for f in os.listdir(BATCH_RUNS_DIR) if isdir(join(BATCH_RUNS_DIR, f)) and len(f.split("_")) > 1], key=lambda x:x.split("_")[-1])
    if not batch_run_dirs:
        raise Exception("Since is_being_processsed returned true there should be at least one batch run")
    return batch_run_dirs[-1]

def get_running_job(config_file):
    config = yaml.load(open(config_file))
    # check that the final deposition step is being processed
    if is_being_processed(config):
    # TODO: check if batched job exists and if so whether last one is running
        run_script = glob.glob(join(get_last_batch_run_dir(), "*{0}".format(BATCH_RUN_SCRIPT_NAME[-3:])))
        if len(run_script) != 1:
            raise Exception("A single run script could not be located: {0}".format(run_script))
        return job_in_queue(run_script[0])
    else:
        return None

def run_batched_deposition(master_config_file, cluster_template_name):
    running_job = get_running_job(master_config_file)
    if not running_job:
        master_config, batch_config_file = setup_batch_config(master_config_file)
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