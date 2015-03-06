import subprocess
import os
from os.path import join, basename, dirname, abspath
from traceback import format_exc
import pmx
import random
import logging
import math
import yaml
import sys
import jinja2
import argparse

VERBOSITY = logging.INFO
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = join(PROJECT_DIR, "templates")

DEBUG = False

GMX_PATH = "/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/bin/"

OUT_STRUCT_FILE = "end.gro"
IN_STRUCT_FILE = "init.gro"

TOPOLOGY_FILE = "topo.top"

TOP_FILE = "templates/{0}.epy".format(TOPOLOGY_FILE)
MDP_FILE = "templates/run.mdp.epy"

GPP_TEMPLATE = "{GMX_PATH}grompp_d -f run.mdp -c {struct} -p topo.top -o md.tpr".format(**{"struct":IN_STRUCT_FILE,
                                                                                                "GMX_PATH":"{GMX_PATH}"})

MPI_ADDITION = "mpirun -np {0} --bind-to-none "

MDRUN = "mdrun_d"
MDRUN_MPI = "mdrun_mpi_d" 

RERUN_FLAG = "-cpi md.cpt -append"


MDRUN_TEMPLATE = "{mpiRun}{GMX_PATH}{mdrun} -pd -s md.tpr -deffnm md -c {struct} {reRunFlag}".format(**{"struct":OUT_STRUCT_FILE,
                                                                                                             "mdrun":       "{mdrun}",
                                                                                                             "GMX_PATH":    "{GMX_PATH}",
                                                                                                             "reRunFlag":   "{reRunFlag}",
                                                                                                             "mpiRun":      "{mpiRun}"}) 

RERUN_SETUP_TEMPLATE = "{GMX_PATH}tpbconv_d -s md.tpr -extend {run_time} -o md.tpr" 


K_B = 0.00831451 #kJ / (mol K)

class Deposition(object):
    
    def __init__(self, runConfigFile, starting_deposition_number=None):
        
        self.runConfig = yaml.load(open(runConfigFile))
        # Convert "*_file" paths in runConfig to be absolute
        recursiveCorrectPaths(self.runConfig, dirname(abspath(runConfigFile)))
        
        # Read starting_deposition_number from YAML file unless provided by command line
        self.moleculeNumber = self.runConfig["starting_deposition_number"] if not starting_deposition_number else starting_deposition_number
        self.rootdir = os.path.abspath(self.runConfig["work_directory"])
        
        # Create the 'work_directory' if it doesn't exist
        if not os.path.exists(self.rootdir):    
            os.makedirs(self.rootdir)
        
        self.rundir = join(self.rootdir, str(self.moleculeNumber))
        
        if self.moleculeNumber == 0:
            configurationPath = self.runConfig["substrate"]["pdb_file"]
        else:
            configurationPath = join(self.rundir, OUT_STRUCT_FILE)
        
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()
       
        self.mixture = None
        self.sampling_mixture = None
        if self.moleculeNumber != 0 :
            self.updateDepositionStep()
            self.countResiduesFromModel()
 
        self.log = "out.log"
        self.err = "out.err"

    def updateDepositionStep(self):
        # Get the list of all the deposition steps
        self.deposition_steps = self.runConfig['deposition_steps']
        # Filter the one that apply to this particular Deposition run
        self.deposition_steps = [ x for x in self.deposition_steps if x["first_sim_id"] <= self.moleculeNumber <= x["last_sim_id"] ]
        # Just to be sure, prevent overlapping deposition step simulation boundaries
        if len(self.deposition_steps) > 1 :
            raise Exception("Overlapping deposition steps definitions. Aborting.")
        elif len(self.deposition_steps) == 0 :
            raise Exception("No deposition step defined for deposition molecule number {0}. Aborting.".format(self.moleculeNumber))

        self.deposition_step = self.deposition_steps[0]
        # Set the sampling mixture
        self.sampling_mixture = self.deposition_step['mixture']
        # Then, update the sampling boundaries with the (maybe new) mixture
        self.setMixtureSamplingBoundaries()
        # Finally, add the potential new residue to the actual system mixture
        if self.mixture :
            # Set the default count to 0 for "new" residues
            new_mixture = self.deposition_step['mixture']
            map(lambda x: x.setdefault("count", 0), new_mixture.values() )
            self.mixture = dict(new_mixture.items() + self.mixture.items()) # Update the mixture dictionnary
        else :
            self.mixture = self.deposition_step['mixture']
        
    def updateModel(self, configurationPath):
        logging.info("Updating model with new configuration")
        #self.model.updateGRO(configurationPath)
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()
        
    def runSystem(self, rerun=False):
        
        reRunFlag = RERUN_FLAG if rerun else ""
        
        if self.runConfig["run_with_mpi"]:
            mpiRun = MPI_ADDITION.format(self.runConfig["max_cores"]) if self.moleculeNumber > self.runConfig["max_cores"] else MPI_ADDITION.format(self.moleculeNumber)
            mdrun = MDRUN_MPI
        else:
            mdrun = MDRUN
            mpiRun = ""
            
        inserts = {"GMX_PATH":   GMX_PATH,
                    "moleculeNumber": self.moleculeNumber,
                    "reRunFlag": reRunFlag,
                    "mpiRun":   mpiRun,
                    "mdrun":    mdrun}
        
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
            errorLog = open(join(self.rootdir, str(self.moleculeNumber), self.err)).readlines()
            for i, l in enumerate(errorLog):
                if "fatal error" in l.lower():
                    logging.error("Error from GROMACS:\n{0}\n".format("".join(errorLog[i:i+7])))
                    break
        
        # now update model after run
        self.updateModel(configurationPath)
    
    def runTPBConf(self, inserts):
        inserts["run_time"] = self.runConfig["run_time"]
        self.run(RERUN_SETUP_TEMPLATE, inserts)
    
    def runGPP(self, inserts):
        
        self.run(GPP_TEMPLATE, inserts)
        
    def runSetup(self):
        self.rundir = join(self.rootdir, str(self.moleculeNumber))
        
        if not os.path.exists(self.rundir):
            os.mkdir(self.rundir)
        
        with open(TOP_FILE) as fh:
            topTemplate = jinja2.Template(fh.read())
            
        with open(MDP_FILE) as fh:
            mdpTemplate = jinja2.Template(fh.read())
        
        with open(join(self.rundir, basename(TOP_FILE)[:-4]),"w") as fh:
            fh.write(topTemplate.render(resMixture=self.mixture, substrate=self.runConfig["substrate"], resnameClusters=cluster(map( lambda x:x.resname, self.model.residues))))
        
        with open(join(self.rundir, basename(MDP_FILE)[:-4]),"w") as fh:
            resList = [res_name for res_name, res in self.mixture.items() if res["count"] > 0]
            fh.write(mdpTemplate.render(resList=resList, substrate=self.runConfig["substrate"], resLength=len(resList), numberOfSteps=int(self.runConfig["run_time"]/self.runConfig["time_step"]), temperature=self.runConfig["temperature"]))
        
        
    
    def run(self, argString, inserts):
        
        argString = argString.format(**inserts)
        argList = argString.split()
        
        logFile = open(join(self.rundir, self.log),"a")
        errFile = open(join(self.rundir, self.err),"a")
        try:
            logging.debug("Running from: '{0}'".format(self.rundir))
            logging.debug("Running command: '{0}'".format(argString))
            subprocess.Popen(argList, cwd=self.rundir, stdout=logFile, stderr=errFile).wait()
        except:
            print "Subprocess terminated with error: \n{0}\n\n{1}".format(argString, format_exc())
        finally: 
            logFile.close()
            errFile.close()
        
        
        
    
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
        logging.debug("Max layer height {0}".format(maxLayerHeight))
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
        if DEBUG :
            highest_atom_height = max([a.x[2] for a in res.atoms])
            logging.debug("Net Z velocity for residue {0}: {1}; Highest Atom Height: {2}".format(res.id, net_z_velocity, highest_atom_height))
        return net_z_velocity >= 0.001 and any([a.x[2] > maxLayerHeight + self.runConfig["escape_tolerance"] for a in res.atoms])
    
    def genInitialVelocitiesLastResidue(self):
        
        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B*self.runConfig["temperature"]/atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            atom.v[2] = -abs(random.gauss(self.runConfig["drift_velocity"], sigma))

    ## this method generates regions between 0-1 that when sampled correspond to particular residues in the mixture
    # e.g. A 1:9 ratio of res1 to res2 is given by: [0.0, 0.1] -> res1, [0.1, 1.0] -> res2
    def setMixtureSamplingBoundaries(self):
        ratioSum = float(sum([v["ratio"] for v in self.sampling_mixture.values()]))
        logging.debug("Ratio sum: " + str(ratioSum))
        movingBoundary = 0.0
        for res in self.sampling_mixture.values():
            res["sample_boundary_min"] = movingBoundary
            movingBoundary += res["ratio"] / ratioSum
            res["sample_boundary_max"] = movingBoundary
        logging.debug("Sampling boundaries: {0}".format(",".join([str((r["res_name"], r["sample_boundary_min"],r["sample_boundary_max"])) for r in self.sampling_mixture.values()])))
 
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
        # Decrease the molecule number
        self.moleculeNumber -= 1
        # Decrease the mixture counts
        self.mixture[residue.resname]["count"] -= 1

    def removeResidue(self, residue):
        self.removeResidueWithID(residue.id)

    def writeInitConfiguration(self):
        updatedPDBPath = join(self.rundir, IN_STRUCT_FILE)
        self.model.write(updatedPDBPath, "{0} deposited molecules".format(self.moleculeNumber), 0)

    def countResiduesFromModel(self):
        for res in self.model.residues:
            # ignore the substrate residue
            if res.resname != self.runConfig["substrate"]["res_name"]: 
                self.mixture[res.resname].setdefault("count", 0)
                self.mixture[res.resname]["count"] += 1
# Why is this needed ?
        # now add in count of zero for any residues not already in the model
        for res in self.mixture.values():
            if not res.has_key("count"):
                res["count"] = 0
# End Why

def recursiveCorrectPaths(node, runConfigFileDir):
    for key, value in node.items():
        if isinstance(value, dict):
            recursiveCorrectPaths(value, runConfigFileDir)
        elif isinstance(value, list):
            for item in value:
                recursiveCorrectPaths(item, runConfigFileDir) 
        elif "file" in key:
            node[key] = join(runConfigFileDir, value)

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

def runDeposition(runConfigFile, starting_deposition_number=None):
    
    deposition = Deposition(runConfigFile, starting_deposition_number=starting_deposition_number)
    
    logging.info("Running deposition with drift velocity of {0:.3f} nm/ps".format(deposition.runConfig["drift_velocity"]))
    while deposition.moleculeNumber < deposition.runConfig["final_deposition_number"]:
        # Increment deposition molecule number
        deposition.moleculeNumber += 1

        # Update the deposition step in case we enter a new deposition phase
        deposition.updateDepositionStep()
        
        # Get the next molecule and insert it into the deposition model with random position, orientation and velocity
        nextMolecule = deposition.getNextMolecule()
        deposition.model.insert_residue(deposition.moleculeNumber, nextMolecule.residues[0], " ")
        deposition.genInitialVelocitiesLastResidue()
     
        # create run directory and run setup make file
        deposition.runSetup()
        
        # write updated model to run directory
        deposition.writeInitConfiguration()
        
        actualMixture = ",".join([" {0}:{1}".format(r["res_name"], r["count"]) for r in deposition.mixture.values()])
        # Do first Run
        logging.info("Running with {0} molecules".format(actualMixture))
        deposition.runSystem()
        
        while not deposition.hasResidueReachedLayer(-1): # -1 means last residue
            if deposition.hasResidueBounced(-1): # -1 means last residue
                logging.warning('It seems like the last inserted molecule has bounced off the surface.')
                deposition.removeResidueWithID(-1) #-1 means last residue
                break
            else:
                logging.info("Rerunning with {0} molecules due to last inserted  molecule not reaching layer".format(actualMixture))
                deposition.runSystem(rerun=True)

        # Iterate over the residues and remove the ones that left the layer
        for residue in deposition.model.residues[1:]: # Dont't try to make sure the substrate is not leaving the layer !
            residue_id = residue.id
            if deposition.hasResidueLeftLayer(residue_id):
                deposition.removeResidueWithID(residue_id)
        
    logging.info("Finished deposition of {0} molecules".format(actualMixture))
    

def parseCommandLine():
    global DEBUG, VERBOSITY
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--debug', dest='debug', action='store_true')
    parser.add_argument('--start', help='{int} Provide a starting deposition number different from the one in the YAML file. Usedd to restart a deposition.')
    args = parser.parse_args()

    if args.debug :
        VERBOSITY = logging.DEBUG
        DEBUG = True
    logging.basicConfig(level=VERBOSITY, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')

    runDeposition(args.input, starting_deposition_number=int(args.start) if args.start else None)
    
if __name__=="__main__":
    parseCommandLine()
