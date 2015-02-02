import subprocess
import os
from os.path import join
from traceback import format_exc
import pmx
import random
import logging

GMX_PATH = "/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/bin/"
OUT_PDB_FILE = "cbp{cbpNumber}-end.gro"
IN_PDB_FILE = "cbp{cbpNumber}-init.gro"

CBP_PDB_000 = "templates/cbp.pdb"
GRAPHENE_PDB = "templates/graphene.pdb"

SETUP_TEMPLATE = "make {cbpNumber}"

GPP_TEMPLATE = "{GMX_PATH}grompp_d -f run.mdp -c {pdbTemplate} -p topo.top -o md.tpr".format(**{"pdbTemplate":IN_PDB_FILE,
                                                                                                "GMX_PATH":"{GMX_PATH}"})

MPI_MDRUN_TEMPLATE = "mpirun -np 8 --bind-to-none {GMX_PATH}mdrun_mpi_d -np 8 -pd -s md.tpr -deffnm md -c {pdbTemplate}".format(**{"pdbTemplate":OUT_PDB_FILE,
                                                                                                                               "GMX_PATH":"{GMX_PATH}"})
 
MDRUN_TEMPLATE = "{GMX_PATH}mdrun_d -pd -s md.tpr -deffnm md -c {pdbTemplate}".format(**{"pdbTemplate":OUT_PDB_FILE,
                                                                                                                               "GMX_PATH":"{GMX_PATH}"}) 


INSERT_DISTANCE = 40.0 #angstrom 
CONTACT_TOLERANCE = 4.0

ROOT_DIR = os.path.abspath("./")


class CBPDeposition(object):
    
    def __init__(self, cbpNumber):
        
        self.cbpNumber = cbpNumber
        
        if cbpNumber == 0:
            pdbPath = join(ROOT_DIR, GRAPHENE_PDB)
        else:
            pdbPath = join(ROOT_DIR, str(self.cbpNumber), OUT_PDB_FILE.format(**{"cbpNumber":self.cbpNumber}))
            
        self.model = pmx.Model(pdbPath)
        
        self.log = "out.log"
        self.err = "out.err"
        
        
    def runSystem(self):
        inserts = {"GMX_PATH": GMX_PATH,
                        "cbpNumber": self.cbpNumber}
        self.run(GPP_TEMPLATE, inserts, join(ROOT_DIR, str(self.cbpNumber)), self.log, self.err)
        self.run(MDRUN_TEMPLATE, inserts, join(ROOT_DIR, str(self.cbpNumber)), self.log, self.err)
    
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
        return any([a.x[2] < maxLayerHeight + CONTACT_TOLERANCE for a in lastRes])

def main(cbpCount, cbpMax):
    
    logging.basicConfig(format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
    
    while cbpCount < cbpMax:
        logging.info("Running with {0} CBP molecules".format(cbpCount))
        
        cBPDeposition = CBPDeposition(cbpCount)
        
        #if not cBPDeposition.hasReachedLayer():
        
        insertHeight = cBPDeposition.getInsertHeight()
        xPos, yPos = cBPDeposition.getRandomPosXY()
        
        cBPModel_000 = pmx.Model(join(ROOT_DIR, CBP_PDB_000))
        
        # set masses to 1.0 to avoid warning
        for atom in cBPModel_000.atoms:
            atom.m = 1.0
            
        cBPModel_000.translate([xPos, yPos, insertHeight])
        cBPModel_000.random_rotation()
        
        
        cBPDeposition.model.insert_residue(cBPDeposition.cbpNumber + 1, cBPModel_000.residues[0], " ")
        cBPDeposition.cbpNumber += 1
        cbpCount += 1
        
        # create run directory
        cBPDeposition.runSetup()
        
        # write updated model to run directory
        updatedPDBPath = join(ROOT_DIR, str(cBPDeposition.cbpNumber), IN_PDB_FILE.format(**{"cbpNumber":cBPDeposition.cbpNumber}))
        cBPDeposition.model.write(updatedPDBPath, "{0} CBP molecules".format(cBPDeposition.cbpNumber), 0)

        cBPDeposition.runSystem()
    
    

main(0, 2)