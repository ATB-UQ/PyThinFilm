import subprocess
import os
from os.path import join
from traceback import format_exc
import pmx
import random
import logging
import math

RUN_WITH_MPI = False

GMX_PATH = "/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/bin/"
OUT_PDB_FILE = "cbp{cbpNumber}-end.gro"
IN_PDB_FILE = "cbp{cbpNumber}-init.gro"

CBP_PDB_000 = "templates/cbp.pdb"
GRAPHENE_PDB = "templates/graphene.pdb"
ITP = "templates/cbp-massmod.itp"

INSERT_DISTANCE = 15.0 #angstrom 
CONTACT_TOLERANCE = 4.0
DRIFT_TIME = 4.0 #ps
DRIFT_VEL = INSERT_DISTANCE*1e-1/DRIFT_TIME


SETUP_TEMPLATE = "make {cbpNumber}"

GPP_TEMPLATE = "{GMX_PATH}grompp_d -f run.mdp -c {pdbTemplate} -p topo.top -o md.tpr".format(**{"pdbTemplate":IN_PDB_FILE,
                                                                                                "GMX_PATH":"{GMX_PATH}"})

MPI_ADDITION = "mpirun -np 16 --bind-to-none "

RERUN_FLAG = "-cpi md.cpt -append"

MDRUN_TEMPLATE = "{mpiRun}{GMX_PATH}mdrun_d -pd -s md.tpr -deffnm md -c {pdbTemplate} {reRunFlag}".format(**{"pdbTemplate":OUT_PDB_FILE,
                                                                                                             "GMX_PATH":    "{GMX_PATH}",
                                                                                                             "reRunFlag":   "{reRunFlag}",
                                                                                                             "mpiRun":      "{mpiRun}"}) 

RERUN_TEMPLATE = "{GMX_PATH}tpbconv_d -s md.tpr -extend {extendTime} -o md.tpr".format(**{"GMX_PATH":"{GMX_PATH}",
                                                                                        "extendTime":DRIFT_TIME}) 

TEMPERATURE = 300 #k

K_B = 0.00831451 #kJ / (mol K)

ROOT_DIR = os.path.abspath("./")


class CBPDeposition(object):
    
    def __init__(self, cbpNumber):
        
        self.cbpNumber = cbpNumber
        
        if cbpNumber == 0:
            pdbPath = join(ROOT_DIR, GRAPHENE_PDB)
        else:
            pdbPath = join(ROOT_DIR, str(self.cbpNumber), OUT_PDB_FILE.format(**{"cbpNumber":self.cbpNumber}))
            
        self.model = pmx.Model(pdbPath)
        self.model.nm2a()
        
        self.log = "out.log"
        self.err = "out.err"
        
        
    def runSystem(self, rerun):
        
        reRunFlag = RERUN_FLAG if rerun else ""
        mpiRun = MPI_ADDITION if RUN_WITH_MPI else ""
        
        inserts = {"GMX_PATH":   GMX_PATH,
                    "cbpNumber": self.cbpNumber,
                    "reRunFlag": reRunFlag,
                    "mpiRun": mpiRun}
        if rerun:
            self.run(RERUN_TEMPLATE, inserts, join(ROOT_DIR, str(self.cbpNumber)), self.log, self.err)
        
        self.run(MDRUN_TEMPLATE, inserts, join(ROOT_DIR, str(self.cbpNumber)), self.log, self.err)
    
    def runGPP(self):
        inserts = {"GMX_PATH":  GMX_PATH,
                   "cbpNumber": self.cbpNumber}
        self.run(GPP_TEMPLATE, inserts, join(ROOT_DIR, str(self.cbpNumber)), self.log, self.err)
        
    def runSetup(self):
        inserts = {"GMX_PATH": GMX_PATH,
                        "cbpNumber": self.cbpNumber}
        self.run(SETUP_TEMPLATE, inserts, ROOT_DIR, self.log, self.err)
    
    def run(self, argString, inserts, workDir, logPath, errPath):
        
        argString = argString.format(**inserts)
        argList = argString.split()
        
        logFile = open(join(workDir, logPath),"a")
        errFile = open(join(workDir, errPath),"a")
        try:
            #print argString
            subprocess.Popen(argList, cwd=workDir, stdout=logFile, stderr=errFile).wait()
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
        
def getMassDict(itpString):
    itpString = itpString.split("[ atoms ]")[1].split("[ bonds ]")[0]
    massDict = {}
    for line in itpString.splitlines():
        if not line or line.startswith(";"):
            continue 
        massDict[line.split()[4]] = float(line.split()[7])
    return massDict

def main(cbpCount, cbpMax):
    
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
    cBPDeposition = CBPDeposition(cbpCount)
    
    logging.info("Running deposition with drift velocity of {0:.3f} nm/ps".format(DRIFT_VEL))
    while cBPDeposition.cbpNumber < cbpMax:
        # increment cbp number
        cBPDeposition.cbpNumber += 1
        
        # create run directory and run setup make file
        cBPDeposition.runSetup()
        insertHeight = cBPDeposition.getInsertHeight()
        xPos, yPos = cBPDeposition.getRandomPosXY()
        
        cBPModel_000 = pmx.Model(join(ROOT_DIR, CBP_PDB_000))
        cBPModel_000.nm2a()
        
        with open(join(ROOT_DIR, ITP),"r") as fh:
            cbpITPString = fh.read()
        massDict = getMassDict(cbpITPString)
        
        # set masses to 1.0 to avoid warning
        for atom in cBPModel_000.atoms:
            atom.m = massDict[atom.name]
            
        cBPModel_000.translate([xPos, yPos, insertHeight])
        cBPModel_000.random_rotation()
        
        
        cBPDeposition.model.insert_residue(cBPDeposition.cbpNumber, cBPModel_000.residues[0], " ")
        
        cBPDeposition.genInitialVelocitiesLastResidue()
     
        # write updated model to run directory
        updatedPDBPath = join(ROOT_DIR, str(cBPDeposition.cbpNumber), IN_PDB_FILE.format(**{"cbpNumber":cBPDeposition.cbpNumber}))
        cBPDeposition.model.write(updatedPDBPath, "{0} CBP molecules".format(cBPDeposition.cbpNumber), 0)
        
        # run grompp
        cBPDeposition.runGPP()
        
        # Do first Run
        logging.info("Running with {0} CBP molecules".format(cBPDeposition.cbpNumber))
        cBPDeposition.runSystem(False)
        while not cBPDeposition.hasReachedLayer():
            logging.info("Rerunning with {0} CBP molecules due to molecule not reaching layer".format(cBPDeposition.cbpNumber))
            cBPDeposition.runSystem(True)
        
    logging.info("Finished deposition of {0} CBP molecules".format(cBPDeposition.cbpNumber))
main(0, 4)
