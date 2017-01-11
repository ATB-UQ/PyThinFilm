import subprocess
import os
from os.path import join, basename, exists
from traceback import format_exc
import pmx
import random
import logging
import math
import yaml
import jinja2
import argparse

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = join(PROJECT_DIR, "templates")

DEBUG = False

if DEBUG:
    logging.basicConfig(level=logging.DEBUG)

OUT_STRUCT_FILE = "end.gro"
IN_STRUCT_FILE = "init.gro"

BASENAME_REMOVE_SUFFIX = lambda path: ".".join( basename(path) .split('.')[0:2])

TOPOLOGY_FILE = "topo.top"

TOP_TEMPLATE = join(TEMPLATE_DIR, "{0}.epy".format(TOPOLOGY_FILE))
TOP_FILE = BASENAME_REMOVE_SUFFIX(TOP_TEMPLATE)
TEMPLATE_ALLOWED_TYPES = ['deposition', 'annealing']

GPP_TEMPLATE = "{{GMX_PATH}}{{grompp}} -f {{MDP_FILE}} -c {struct} -p topo.top -o md.tpr".format(struct=IN_STRUCT_FILE)

MPI_ADDITION = "mpirun -np {0} --bind-to-none "

GROMPP = "grompp_d"
TPBCONV = "tpbconv_d"
MDRUN = "mdrun_d"
MDRUN_MPI = "mdrun_mpi_d" 

RERUN_FLAG = "-cpi md.cpt -append"


MDRUN_TEMPLATE = "{{mpiRun}}{{GMX_PATH}}{{mdrun}} -pd -s md.tpr -deffnm md -c {struct} {{reRunFlag}}".format(struct=OUT_STRUCT_FILE)

RERUN_SETUP_TEMPLATE = "{GMX_PATH}{tpbconv} -s md.tpr -extend {run_time} -o md.tpr" 


K_B = 0.00831451 #kJ / (mol K)

DEFAULT_PARAMETERS = {"lincs_order": 4,
              "lincs_iterations":1}

class Deposition(object):

    def __init__(self, runConfigFile, starting_deposition_number=None):
        self.runConfig = yaml.load(open(runConfigFile))
        # Convert template paths in runConfig to be absolute
        recursiveCorrectPaths(self.runConfig, PROJECT_DIR)

        # Read starting_deposition_number from YAML file unless provided by command line
        self.run_ID = self.runConfig["starting_deposition_number"] if not starting_deposition_number else starting_deposition_number
        self.rootdir = os.path.abspath(self.runConfig["work_directory"])
        self.disambiguate_run_config()

        # Create the 'work_directory' if it doesn't exist
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)

        # By default, try restarting form the current branch
        self.rundir = self.get_rundir()
        # Otherwise restart from the master '' branch
        if not exists(self.rundir):
            self.rundir = join(self.rootdir, str(self.run_ID))
        print self.rundir

        if self.run_ID == 0:
            configurationPath = self.runConfig["substrate"]["pdb_file"]
        else:
            configurationPath = join(self.rundir, OUT_STRUCT_FILE)

        self.model = pmx.Model(configurationPath)
        self.model.nm2a()

        self.mixture = {}
        self.sampling_mixture = {}

        self.initDepositionSteps()

        self.log = "out.log"
        self.err = "out.err"

    def molecule_number(self):
        return len(self.model.residues)

    def get_rundir(self):
        suffix = str(self.runConfig['development']['branch']) if ("development" in self.runConfig) else ""
        return join(self.rootdir, str(self.run_ID) + suffix)

    def initDepositionSteps(self):
        # Get the list of all the deposition steps
        self.deposition_steps = self.runConfig['deposition_steps']
        self.last_run_ID = self.runConfig["final_deposition_number"]
        self.initMixtureAndResidueCounts()

    def disambiguate_run_config(self):
        self.gmx_path = self.runConfig['gmx_path'] if self.runConfig['gmx_path'] else ""

    def updateDepositionStep(self):
        # Filter the one that apply to this particular Deposition run
        self.current_deposition_steps = [ x for x in self.deposition_steps if x["first_sim_id"] <= self.run_ID <= x["last_sim_id"] ]
        # Just to be sure, prevent overlapping deposition step simulation boundaries
        if len(self.current_deposition_steps) > 1 :
            raise Exception("Overlapping deposition steps definitions. Aborting.")
        elif len(self.current_deposition_steps) == 0 :
            raise Exception("No deposition step defined for deposition molecule number {0}. Aborting.".format(self.run_ID))

        self.deposition_step = self.current_deposition_steps[0]
        # Update the MDP template
        self.template_type = self.deposition_step['template']['type']
        if not self.template_type in TEMPLATE_ALLOWED_TYPES :
            raise Exception('Unknown template type: {template_type}. Allowed types are: {allowed_types}'.format(template_type=self.template_type, allowed_types=TEMPLATE_ALLOWED_TYPES) )
        self.mdp_template_file = self.deposition_step['template']['file']
        self.mdp_file = BASENAME_REMOVE_SUFFIX(self.mdp_template_file)

        # For depositon steps only, update the sampling mixture
        if self.isDepositionRun():
            # Set the sampling mixture to the next current step mixture
            self.sampling_mixture = self.deposition_step['mixture']
            # Then, update the sampling boundaries with the (maybe new) mixture
            self.setMixtureSamplingBoundaries()

    def isDepositionRun(self):
        return self.template_type == 'deposition'
    def isAnnealingRun(self):
        return self.template_type == 'annealing'

    def runParameters(self):
        parameters_dict =  { \
                'run_ID': "{0}/{1}".format(self.run_ID, self.last_run_ID),
                'temperature':    "{0} K".format(self.deposition_step["temperature"]),
                'run_time':       "{0} ps".format(self.deposition_step["run_time"]),
        }
        if self.isDepositionRun():
            parameters_dict['drift_velocity'] = "{0} nm/ps".format(self.runConfig["drift_velocity"])
        return parameters_dict

    def updateModel(self, configurationPath):
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()

    def runSystem(self, rerun=False):

        reRunFlag = RERUN_FLAG if rerun else ""

        if self.runConfig["run_with_mpi"]:
            mpiRun = MPI_ADDITION.format(self.runConfig["max_cores"]) if self.molecule_number() > self.runConfig["max_cores"] else MPI_ADDITION.format(self.molecule_number())
            mdrun = self.runConfig["mdrun_mpi"] \
                        if "mdrun_mpi" in self.runConfig else MDRUN_MPI
        else:
            mdrun = self.runConfig["mdrun"] \
                        if "mdrun" in self.runConfig else MDRUN
            mpiRun = ""
        grompp = self.runConfig["grompp"] \
                    if "grompp" in self.runConfig else GROMPP

        inserts = {"GMX_PATH":   self.gmx_path,
                   "MDP_FILE": self.mdp_file,
                   "run_ID": self.run_ID,
                   "reRunFlag": reRunFlag,
                   "mpiRun":   mpiRun,
                   "mdrun":    mdrun, 
                   "grompp": grompp,
                   }

        if rerun:
            # run rerun setup script
            self.runTPBConf(inserts)
        else:
            # run grompp 
            self.runGPP(inserts)

        # run the md
        self.run(MDRUN_TEMPLATE, inserts)

        configurationPath = join(self.rundir, OUT_STRUCT_FILE) 

        if not os.path.exists(configurationPath):
            logging.error("MD run did not produce expected output file")
            errorLog = open(join(self.rootdir, str(self.run_ID), self.err)).readlines()
            for i, l in enumerate(errorLog):
                if "fatal error" in l.lower():
                    logging.error("Error from GROMACS:\n{0}\n".format("".join(errorLog[i:i+7])))
                    break

        # now update model after run
        self.updateModel(configurationPath)

    def runTPBConf(self, inserts):
        inserts["run_time"] = self.deposition_step["run_time"]
        inserts["tpbconv"] = self.runConfig["tpbconv"] \
                    if "tpbconv" in self.runConfig else TPBCONV
        self.run(RERUN_SETUP_TEMPLATE, inserts)

    def runGPP(self, inserts):

        self.run(GPP_TEMPLATE, inserts)

    def runSetup(self):
        self.rundir = self.get_rundir()

        if not os.path.exists(self.rundir):
            os.mkdir(self.rundir)

        with open(TOP_TEMPLATE) as fh:
            topTemplate = jinja2.Template(fh.read())

        with open(self.mdp_template_file) as fh:
            mdpTemplate = jinja2.Template(fh.read())

        with open(join(self.rundir, TOP_FILE), "w") as fh:
            fh.write(topTemplate.render(resMixture=self.mixture, substrate=self.runConfig["substrate"], resnameClusters=cluster(map( lambda x:x.resname, self.model.residues))))

        with open(join(self.rundir, self.mdp_file),"w") as fh:
            resList = [res_name for res_name, res in self.mixture.items() if res["count"] > 0]
            fh.write(mdpTemplate.render(resList=resList, 
                    substrate=self.runConfig["substrate"], 
                    resLength=len(resList), 
                    timeStep=self.runConfig["time_step"], 
                    numberOfSteps=int(self.deposition_step["run_time"]/self.runConfig["time_step"]), 
                    temperature=self.deposition_step["temperature"],
                    lincs_order=self.deposition_step["lincs_order"] if "lincs_order" in self.deposition_step else DEFAULT_PARAMETERS["lincs_order"],
                    lincs_iterations= self.deposition_step["lincs_iterations"] if "lincs_iterations" in self.deposition_step else DEFAULT_PARAMETERS["lincs_iterations"],
                    )
            )

    def run(self, argString, inserts):

        argString = argString.format(**inserts)
        argList = argString.split()

        returncode = 0
        logFile = open(join(self.rundir, self.log),"a")
        errFile = open(join(self.rundir, self.err),"a")
        try:
            logging.debug("    Running from: '{0}'".format(self.rundir))
            logging.debug("    Running command: '{0}'".format(argString))
            returncode = subprocess.Popen(argList, cwd=self.rundir, stdout=logFile, stderr=errFile, env=os.environ).wait()
        except:
            logging.error("Subprocess terminated with error: \n{0}\n\n{1}".format(argString, format_exc()))
            raise
        finally: 
            logFile.close()
            errFile.close()
        if 0 < returncode:
            msg = "Subprocess terminated with nonzero exit code: \n{0}\n\n{1}".format(argString, format_exc())
            logging.error(msg)
            raise Exception(msg)

    def getInsertHeight(self):
        return  self.maxZHeight() + self.runConfig["insert_distance"]

    def maxZHeight(self):
        return max([a.x[2] for a in self.model.atoms])

    def getRandomPosXY(self):
        x, y = map(lambda x:[v*10 for v in x], self.model.box[:2])
        # take the 0th and 1st element from x and y respectively as these must be rectangular boxes
        return random.uniform(0.0, x[0]), random.uniform(0.0, y[1])

    def maxLayerHeight(self, excluded_res):
        maxLayerHeight = max([a.x[2] for a in self.model.atoms if a not in excluded_res.atoms])
        logging.debug("    Max layer height {0}".format(maxLayerHeight))
        return maxLayerHeight

    def hasResidueReachedLayer(self, residue_ID):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        maxLayerHeight = self.maxLayerHeight(res)
        return any([a.x[2] < maxLayerHeight + self.runConfig["contact_tolerance"] for a in res.atoms])

    # A molecule has bounced if any of its atoms is above its insertion height
    # NB: Relies on the definition of self.insertionHeight
    def hasResidueBounced(self, residue_ID):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        return any([a.x[2] > self.insertionHeight for a in res.atoms])

    # A molecule as left layer if it is above a certain height (above the layer's mean height ??)  with a net z velocity
    def hasResidueLeftLayer(self, residue_ID):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        maxLayerHeight = self.maxLayerHeight(res)
        #  Should it be mass weighted ?
        net_z_velocity = sum( map( lambda x: x.v[2], res.atoms) ) / len(res.atoms)
        highest_atom_height = max([a.x[2] for a in res.atoms])
        logging.debug("    Net Z velocity for residue {0}: {1}; Highest Atom Height: {2}".format(res.id, net_z_velocity, highest_atom_height))
        return net_z_velocity > 0. and any([a.x[2] > maxLayerHeight + self.runConfig["escape_tolerance"] for a in res.atoms])

    def genInitialVelocitiesLastResidue(self):

        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B*self.deposition_step["temperature"]/atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            atom.v[2] = -abs(random.gauss(self.runConfig["drift_velocity"], sigma))

    ## this method generates regions between 0-1 that when sampled correspond to particular residues in the mixture
    # e.g. A 1:9 ratio of res1 to res2 is given by: [0.0, 0.1] -> res1, [0.1, 1.0] -> res2
    def setMixtureSamplingBoundaries(self):
        ratioSum = float(sum([v["ratio"] for v in self.sampling_mixture.values()]))
        movingBoundary = 0.0
        for res in self.sampling_mixture.values():
            res["sample_boundary_min"] = movingBoundary
            movingBoundary += res["ratio"] / ratioSum
            res["sample_boundary_max"] = movingBoundary

    def sampleMixture(self):
        if len(self.sampling_mixture) == 1:
            return self.sampling_mixture.values()[0]
        else:
            randomNumber = random.uniform(0.0,1.0)
            for res in self.sampling_mixture.values():
                if res["sample_boundary_min"] <= randomNumber <= res["sample_boundary_max"]:
                    return res


    def getNextMolecule(self):
        nextMol = self.sampleMixture()
        nextMolecule = pmx.Model(nextMol["pdb_file"])
        nextMolecule.nm2a()

        self.mixture[nextMol["res_name"]].setdefault("count", 0)
        self.mixture[nextMol["res_name"]]["count"] += 1
        logging.info("    Inserting molecule: {0}".format(nextMol["res_name"]) )

        with open(nextMol["itp_file"]) as fh:
            cbpITPString = fh.read()
        massDict = getMassDict(cbpITPString)

        # set masses to 1.0 to avoid warning
        for atom in nextMolecule.atoms:
            atom.m = massDict[atom.name]

        insertHeight = self.getInsertHeight()
        xPos, yPos = self.getRandomPosXY()

        nextMolecule.translate([xPos, yPos, insertHeight])
        nextMolecule.random_rotation()

        # Define its insertion height as the highest Z value amongst its atoms
        self.insertionHeight = max([ a.x[2] for a in nextMolecule.atoms])

        return nextMolecule

    def removeResidueWithID(self, residue_ID):
        residue = self.model.residues[residue_ID]
        self.model.remove_residue(residue)
        logging.debug("Removing residue: {0}".format(residue.resname))
        # Decrease the run ID by one
        self.run_ID -= 1
        # Decrease the mixture counts
        self.mixture[residue.resname]["count"] -= 1

    def removeResidue(self, residue):
        self.removeResidueWithID(residue.id)

    def writeInitConfiguration(self):
        updatedPDBPath = join(self.rundir, IN_STRUCT_FILE)
        self.model.write(updatedPDBPath, "{0} deposited molecules".format(self.run_ID), 0)

    def initMixtureAndResidueCounts(self):
        # Set the count to zero for all the residues in the different mixtures from the different deposition steps
        for mixture in [ x['mixture'] for x in self.deposition_steps if 'mixture' in x ] :
            for res in mixture.values():
                resname = res["res_name"]
                self.mixture[resname] = res
                self.mixture[resname]["count"]= 0

        # Then, update count for residues already in the model
        for res in self.model.residues:
            # ignore the substrate residue
            if res.resname != self.runConfig["substrate"]["res_name"]: 
                self.mixture[res.resname]["count"] += 1
        # Log the mixture
        logging.debug('Initial mixture is: {mixture}'.format(mixture=self.mixture))

def recursiveCorrectPaths(node, root_dir):
    for key, value in node.items():
        if isinstance(value, dict):
            recursiveCorrectPaths(value, root_dir)
        elif isinstance(value, list):
            for item in value:
                recursiveCorrectPaths(item, root_dir) 
        elif "file" in key:
            node[key] = join(root_dir, value)

def getMassDict(itpString):
    itpString = itpString.split("[ atoms ]")[1].split("[ bonds ]")[0]
    massDict = {}
    for line in itpString.splitlines():
        if not line or line.startswith(";"):
            continue 
        massDict[line.split()[4]] = float(line.split()[7])
    return massDict

def cluster(resnameList):
    clusterList = []
    i = 0
    current_resname = resnameList[0]
    if not resnameList:
        return ""
    while current_resname != '':
        count = 1
        try:
            next_resname = resnameList[i+1]
        except IndexError:
            next_resname = ""
        while next_resname == current_resname :
            count +=1
            try:
                next_resname = resnameList[i+count]
            except IndexError:
                next_resname = ""
        clusterList.append("{0} {1}".format(current_resname, count))
        i += count
        try:
            current_resname = resnameList[i]
        except IndexError:
            current_resname = ""
    return clusterList

def runDeposition(runConfigFile, starting_deposition_number=None, remove_bounce=False, remove_leaving_layer=False, debug=DEBUG):
    if debug:
        verbosity = logging.DEBUG
        format_log = '%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)'
    else:
        verbosity = logging.INFO
        format_log = '%(asctime)s - [%(levelname)s] - %(message)s'
    logging.basicConfig(level=verbosity, format=format_log, datefmt='%d-%m-%Y %H:%M:%S')
    deposition = Deposition(runConfigFile, starting_deposition_number=starting_deposition_number)

    while deposition.run_ID < deposition.last_run_ID:
        # Increment run ID
        deposition.run_ID += 1

        # Update the deposition step in case we enter a new deposition phase
        deposition.updateDepositionStep()

        # Depending on run type ...
        if deposition.isDepositionRun():
            logging.info("[DEPOSITION] Running deposition run with parameters: {parameters_dict}".format(parameters_dict=deposition.runParameters()))
            # Get the next molecule and insert it into the deposition model with random position, orientation and velocity
            nextMolecule = deposition.getNextMolecule()
            last_residue = len(deposition.model.residues)
            deposition.model.insert_residue(last_residue, nextMolecule.residues[0], " ")
            deposition.genInitialVelocitiesLastResidue()
        elif deposition.isAnnealingRun():
            logging.info("[ANNEALING] Running annealing run with parameters: {parameters_dict}".format(parameters_dict=deposition.runParameters()))

        # create run directory and run setup make file
        deposition.runSetup()

        # Write updated model to run directory
        deposition.writeInitConfiguration()

        actualMixture = ",".join([" {0}:{1}".format(r["res_name"], r["count"]) for r in deposition.mixture.values()])
        # Do first Run
        logging.info("    Current mixture is: {0}".format(actualMixture))
        deposition.runSystem()

        if deposition.isDepositionRun():
            while not deposition.hasResidueReachedLayer(-1): # -1 means last residue
                if remove_bounce and deposition.hasResidueBounced(-1): # -1 means last residue
                    logging.warning('It seems like the last inserted molecule has bounced off the surface.')
                    deposition.removeResidueWithID(-1) #-1 means last residue
                    break
                else:
                    logging.info("    Rerunning with same parameters ({parameters_dict}) due to last inserted  molecule not reaching layer".format(parameters_dict=deposition.runParameters()))
                    deposition.runSystem(rerun=True)

            if remove_leaving_layer :
                # Iterate over the residues and remove the ones that left the layer
                for residue in deposition.model.residues[1:]: # Dont't try to make sure the substrate is not leaving the layer !
                    residue_id = residue.id
                    if deposition.hasResidueLeftLayer(residue_id):
                        deposition.removeResidueWithID(residue_id)

    logging.info("Finished deposition of {0} molecules".format(deposition.last_run_ID))


def parseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--debug', dest='debug', action='store_true')
    parser.add_argument('--start', help='{int} Provide a starting deposition number different from the one in the YAML file. Used to restart a deposition. This number corresponds to the last successful deposition. Use 0 to start from scratch.')
    parser.add_argument('--remove-mol-bouncing', dest='remove_bounce',        action='store_true')
    parser.add_argument('--remove-mol-leaving',  dest='remove_leaving_layer', action='store_true')
    args = parser.parse_args()

    runDeposition(args.input, starting_deposition_number=int(args.start) if args.start else None, remove_bounce=args.remove_bounce, remove_leaving_layer=args.remove_leaving_layer, debug=args.debug)

if __name__=="__main__":
    parseCommandLine()
