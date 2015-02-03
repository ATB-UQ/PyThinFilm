import subprocess
import os
from os.path import join
from traceback import format_exc
import pmx
import random
import logging
import math
import shutil

RUN_WITH_MPI = True
VERBOCITY = logging.INFO
VERBOCITY = logging.DEBUG

MAX_CORES = 8

INSERT_DISTANCE = 20.0 #angstrom 
CONTACT_TOLERANCE = 4.0
DRIFT_TIME = 4.0 #ps

DRIFT_VEL = INSERT_DISTANCE*1e-1/DRIFT_TIME
DRIFT_VEL_REALISTIC = 0.08
DRIFT_VEL = DRIFT_VEL_REALISTIC


PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

GMX_PATH = "/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/bin/"
OUT_STRUCT_FILE = "cbp{moleculeNumber}-end.gro"
IN_STRUCT_FILE = "cbp{moleculeNumber}-init.gro"

TEMPLATE_DIR = join(PROJECT_DIR, "templates")

CBP_PDB_000 = join(TEMPLATE_DIR, "cbp.pdb")

GRAPHENE_PDB = join(TEMPLATE_DIR, "graphene.pdb")
ITP = join(TEMPLATE_DIR, "cbp-massmod.itp")

SETUP_TEMPLATE = "make {moleculeNumber} template_dir={templatePath}"

GPP_TEMPLATE = "{GMX_PATH}grompp_d -f run.mdp -c {pdbTemplate} -p topo.top -o md.tpr".format(**{"pdbTemplate":IN_STRUCT_FILE,
                                                                                                "GMX_PATH":"{GMX_PATH}"})

MPI_ADDITION = "mpirun -np {0} --bind-to-none "

MDRUN = "mdrun_mpi_d" if RUN_WITH_MPI else "mdrun_d"

RERUN_FLAG = "-cpi md.cpt -append"


MDRUN_TEMPLATE = "{mpiRun}{GMX_PATH}{mdrun} -pd -s md.tpr -deffnm md -c {pdbTemplate} {reRunFlag}".format(**{"pdbTemplate":OUT_STRUCT_FILE,
                                                                                                             "mdrun":MDRUN,
                                                                                                             "GMX_PATH":    "{GMX_PATH}",
                                                                                                             "reRunFlag":   "{reRunFlag}",
                                                                                                             "mpiRun":      "{mpiRun}"}) 

RERUN_SETUP_TEMPLATE = "{GMX_PATH}tpbconv_d -s md.tpr -extend {extendTime} -o md.tpr".format(**{"GMX_PATH":"{GMX_PATH}",
                                                                                        "extendTime":DRIFT_TIME}) 

TEMPERATURE = 300 #k

K_B = 0.00831451 #kJ / (mol K)

MIXTURE = {CBP_PDB_000:21,
           }

class Deposition(object):
    
    def __init__(self, moleculeNumber, rootdir):
        
        logging.basicConfig(level=VERBOCITY, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
        
        self.moleculeNumber = moleculeNumber
        self.rootdir = os.path.abspath(rootdir)
        
        if not os.path.exists(self.rootdir):    
            os.makedirs(self.rootdir)
        
        # copy makefile to rootdir
        makefilePath = join(rootdir, "Makefile")
        if not os.path.exists(makefilePath):
            shutil.copy(join(PROJECT_DIR, "Makefile"), makefilePath)
        
        if self.moleculeNumber == 0:
            configurationPath = join(PROJECT_DIR, GRAPHENE_PDB)
        else:
            configurationPath = join(self.rootdir, str(self.moleculeNumber), OUT_STRUCT_FILE.format(**{"moleculeNumber":self.moleculeNumber}))
        
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()
        
        self.log = "out.log"
        self.err = "out.err"
        
        
    def updateModel(self, configurationPath):
        logging.info("Updating model with new configuration")
        #self.model.updateGRO(configurationPath)
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()
        
    def runSystem(self, rerun=False):
        
        reRunFlag = RERUN_FLAG if rerun else ""
        
        if RUN_WITH_MPI:
            mpiRun = MPI_ADDITION.format(MAX_CORES) if self.moleculeNumber > MAX_CORES else MPI_ADDITION.format(self.moleculeNumber) 
        else:
            mpiRun = ""
            
        inserts = {"GMX_PATH":   GMX_PATH,
                    "moleculeNumber": self.moleculeNumber,
                    "reRunFlag": reRunFlag,
                    "mpiRun": mpiRun}
        
        if rerun:
            # run rerun setup script
            self.runTPBConf(inserts)
        else:
            # run grompp 
            self.runGPP(inserts)
        
        # run the md
        self.run(MDRUN_TEMPLATE, inserts)
        
        configurationPath = join(self.rootdir, str(self.moleculeNumber), OUT_STRUCT_FILE.format(**{"moleculeNumber":self.moleculeNumber})) 
        
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
        self.run(RERUN_SETUP_TEMPLATE, inserts)
    
    def runGPP(self, inserts):
        
        self.run(GPP_TEMPLATE, inserts)
        
    def runSetup(self):
        inserts = {"GMX_PATH": GMX_PATH,
                    "moleculeNumber": self.moleculeNumber,
                    "templatePath": TEMPLATE_DIR}
        
        self.run(SETUP_TEMPLATE, inserts, workdir=self.rootdir)
    
    def run(self, argString, inserts, workdir=None):
        
        argString = argString.format(**inserts)
        argList = argString.split()
        if not workdir:
            workdir = join(self.rootdir, str(self.moleculeNumber))
        
        logFile = open(join(workdir, self.log),"a")
        errFile = open(join(workdir, self.err),"a")
        try:
            logging.debug("Running from: '{0}'".format(workdir))
            logging.debug("Running command: '{0}'".format(argString))
            subprocess.Popen(argList, cwd=workdir, stdout=logFile, stderr=errFile).wait()
        except:
            print "Subprocess terminated with error: \n{0}\n\n{1}".format(argString, format_exc())
        finally: 
            logFile.close()
            errFile.close()
        
        
        
    
    def getInsertHeight(self):
        return  self.maxZHeight() + INSERT_DISTANCE
    
    def maxZHeight(self):
        return max([a.x[2] for a in self.model.atoms])
        
    def getRandomPosXY(self):
        x, y = map(lambda x:[v*10 for v in x], self.model.box[:2])
        # take the 0th and 1st element from x and y respectively as these must be rectangular boxes
        return random.uniform(0.0, x[0]), random.uniform(0.0, y[1])

    def hasReachedLayer(self):
        lastResID = self.model.residues[-1].id
        lastRes = self.model.residue(lastResID)
        maxLayerHeight = max([a.x[2] for a in self.model.atoms if a not in lastRes.atoms])
        logging.debug("Max layer height {0}".format(maxLayerHeight))
        return any([a.x[2] < maxLayerHeight + CONTACT_TOLERANCE for a in lastRes.atoms])
    
    def genInitialVelocitiesLastResidue(self):
        
        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B*TEMPERATURE/atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            atom.v[2] = -abs(random.gauss(DRIFT_VEL, sigma))
            
        
    def sampleMixture(self):
        if len(MIXTURE) == 1:
            return join(PROJECT_DIR, MIXTURE.keys()[0])
        
                    
    def getNextMolecule(self):
        nextMolPath = self.sampleMixture()
        nextMolecule = pmx.Model(nextMolPath)
        nextMolecule.nm2a()
        
        with open(join(PROJECT_DIR, ITP),"r") as fh:
            cbpITPString = fh.read()
        massDict = getMassDict(cbpITPString)
        
        # set masses to 1.0 to avoid warning
        for atom in nextMolecule.atoms:
            atom.m = massDict[atom.name]
        
        insertHeight = self.getInsertHeight()
        xPos, yPos = self.getRandomPosXY()
            
        nextMolecule.translate([xPos, yPos, insertHeight])
        nextMolecule.random_rotation()
        return nextMolecule
        
    def writeInitConfiguration(self):
        updatedPDBPath = join(self.rootdir, str(self.moleculeNumber), IN_STRUCT_FILE.format(**{"moleculeNumber":self.moleculeNumber}))
        self.model.write(updatedPDBPath, "{0} CBP molecules".format(self.moleculeNumber), 0)
        
def getMassDict(itpString):
    itpString = itpString.split("[ atoms ]")[1].split("[ bonds ]")[0]
    massDict = {}
    for line in itpString.splitlines():
        if not line or line.startswith(";"):
            continue 
        massDict[line.split()[4]] = float(line.split()[7])
    return massDict

def runDeposition(cbpCount, cbpMax, rootdir):
    
    deposition = Deposition(cbpCount, rootdir)
    
    logging.info("Running deposition with drift velocity of {0:.3f} nm/ps".format(DRIFT_VEL))
    while deposition.moleculeNumber < cbpMax:
        # increment cbp number
        deposition.moleculeNumber += 1
        
        # create run directory and run setup make file
        deposition.runSetup()
        
        # get the next molecule and insert it into the deposition model with random position, orientation and velocity
        nextMolecule = deposition.getNextMolecule()
        deposition.model.insert_residue(deposition.moleculeNumber, nextMolecule.residues[0], " ")
        deposition.genInitialVelocitiesLastResidue()
     
        # write updated model to run directory
        deposition.writeInitConfiguration()
        
        # Do first Run
        logging.info("Running with {0} CBP molecules".format(deposition.moleculeNumber))
        deposition.runSystem()
        
        while not deposition.hasReachedLayer():
            logging.info("Rerunning with {0} CBP molecules due to molecule not reaching layer".format(deposition.moleculeNumber))
            deposition.runSystem(rerun=True)
        
    logging.info("Finished deposition of {0} CBP molecules".format(deposition.moleculeNumber))
    
    
runDeposition(45, 301, os.path.abspath("./dep60"))
#runDeposition(108, 111, os.path.abspath("./firstRun110Mols"))
