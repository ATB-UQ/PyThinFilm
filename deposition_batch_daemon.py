import yaml
import argparse
import os
from os.path import dirname, exists, abspath, basename, join, isdir 
import subprocess
import logging
import glob
from copy import deepcopy

CLUSTER_TEMPLATE_NAME = "{name}_run_script.sh.template"

LOG_FILE = "stdout.log"
BATCH_RUNS_DIR = "batch_runs"
BATCH_RUN_DIR_NAME = "from_{start}_to_{end}"
BATCH_RUN_SCRIPT_NAME = "deposition_batch_{start}_{end}.sh"
BATCH_RUN_CONFIG_NAME = "config.yml"
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

def setup_batch_config(config_file, last_run_dir):
    config = yaml.load(open(config_file))
    last_run = int(basename(last_run_dir))
    last_deposition_steps = get_last_deposition_step(config)

    last_deposition_steps["first_sim_id"] = last_run + 1
    last_deposition_steps["last_sim_id"] = last_run + last_deposition_steps["n_per_batch"]
    last_deposition_steps["description"] = "Batch run from {0} to {1} of: {2}".format(last_deposition_steps["first_sim_id"], last_deposition_steps["last_sim_id"], last_deposition_steps["description"])

    batch_config = deepcopy(config)
    batch_config["deposition_steps"] = [last_deposition_steps]

    batch_run_dir = get_batch_run_dir(last_deposition_steps["first_sim_id"], last_deposition_steps["last_sim_id"])
    if exists(batch_run_dir):
        raise Exception("Batch run already exists: {0}".format(batch_run_dir))
    os.makedirs(batch_run_dir)
    batch_config_file = join(batch_run_dir, BATCH_RUN_CONFIG_NAME)
    with open(batch_config_file, "w") as fh:
        yaml.dump(batch_config, fh)
    return batch_config_file, last_deposition_steps["first_sim_id"], last_deposition_steps["last_sim_id"]

def get_batch_run_dir(start, end):
    return join(BATCH_RUNS_DIR, BATCH_RUN_DIR_NAME.format(start=start, end=end))

def populate_run_script_template(template_name, batch_config_file, start, end):
    cluster_template = open(join("templates", CLUSTER_TEMPLATE_NAME.format(name=template_name))).read()
    batch_run_dir = get_batch_run_dir(start, end)
    log_file = join(batch_run_dir, LOG_FILE)
    run_script = cluster_template.format(config_file=batch_config_file, start_from=start, stdout_file=log_file)
    run_script_file = join(batch_run_dir, BATCH_RUN_SCRIPT_NAME.format(start=start, end=end))
    with open(run_script_file, "w") as fh:
        fh.write(run_script)
    return run_script_file

def submit_to_queue(run_script):
    run(["qsub", run_script], cwd_opt=dirname(abspath(run_script)))

def job_in_queue(run_script):
    return
    qstat = run_with_shell_stdout_piped(PBS_GET_ALL_MY_JOB_NAMES)
    if basename(run_script) in qstat:
        return basename(run_script)


def get_last_deposition_step(config):
    return sorted(config["deposition_steps"], key=lambda x:x["first_sim_id"])[-1]

def is_being_processed(config):
    last_run_dir = get_last_run_dir(config)
    logging.info("Last run dir: {0}".format(last_run_dir))
    last_deposition_steps = get_last_deposition_step(config)
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

def get_last_run_dir(config):
    run_dirs = [join(config["work_directory"], f) for f in os.listdir(config["work_directory"]) if isdir(join(config["work_directory"], f))]
    return sorted(run_dirs)[-1] if run_dirs else None

def get_last_batch_run_dir():
    batch_run_dirs = sorted([join(BATCH_RUNS_DIR, f) for f in os.listdir(BATCH_RUNS_DIR) if isdir(join(BATCH_RUNS_DIR, f)) and len(f.split("_")) > 1], key=lambda x:x.split("_")[-1])
    if not batch_run_dirs:
        raise Exception("Since is_being_processsed returned true there should be at least one batch run")
    return batch_run_dirs[-1]

def get_current_status(config_file):
    config = yaml.load(open(config_file))
    # check that the final deposition is being processed
    last_run_dir = get_last_run_dir(config)
    logging.info("Last run deposition: {0}".format(last_run_dir))
    if is_being_processed(config):
        # check if batched job exists and if so whether last one is running
        run_script = glob.glob(join(get_last_batch_run_dir(), "*{0}".format(BATCH_RUN_SCRIPT_NAME[-3:])))
        if len(run_script) != 1:
            raise Exception("A single run script could not be located: {0}".format(run_script))
        return job_in_queue(run_script[0]), last_run_dir
    else:
        return None, last_run_dir

def run_batched_deposition(master_config_file, cluster_template_name):
    running_job, last_run_dir = get_current_status(master_config_file)
    if not running_job:
        batch_config_file, start, end = setup_batch_config(master_config_file, last_run_dir)
        run_script_file = populate_run_script_template(cluster_template_name, batch_config_file, start, end)
        submit_to_queue(run_script_file)
    else:
        logging.info("A batch job is already running: {0} (currently processing {1})".format(running_job, basename(last_run_dir)))

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--master_config')
    parser.add_argument('-c', '--cluster_template')
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args = parse_command_line()
    run_batched_deposition(args.master_config, args.cluster_template)