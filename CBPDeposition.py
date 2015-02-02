import subprocess
import os
from os.path import join
from traceback import format_exc
import pmx
import random

GMX_PATH = "/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/bin/"
OUT_PDB_FILE = "cbp{cbpNumber}-end.pdb"
IN_PDB_FILE = "cbp{cbpNumber}-init.pdb"

CBP_PDB_000 = "templates/cbp.pdb"

SETUP_TEMPLATE = "make {cbpNumber}"
GPP_TEMPLATE = "{GMX_PATH}grompp_d -f run.mdp -c {pdbTemplate} -p topo.top -o md.tpr".format(**{"pdbTemplate":IN_PDB_FILE,
                                                                                                "GMX_PATH":"{GMX_PATH}"})

MDRUN_TEMPLATE = "mpirun -np 8 --bind-to-none {GMX_PATH}mdrun_mpi_d -np 8 -pd -s md.tpr -deffnm md -c {pdbTemplate}".format(**{"pdbTemplate":OUT_PDB_FILE,
                                                                                                                               "GMX_PATH":"{GMX_PATH}"}) 

INSERT_DISTANCE = 40.0 #angstrom 
CONTACT_TOLERANCE = 4.0

ROOT_DIR = os.path.abspath("./")


class CBPDeposition(object):
    
    def __init__(self, cbpNumber):
        
        self.cbpNumber = cbpNumber
        
        pdbPath = join(ROOT_DIR, self.cbpNumber, OUT_PDB_FILE.format(**{"cbpNumber":self.cbpNumber}))
        self.model = pmx.Model(pdbPath)
        
        
    def runSystem(self):
        inserts = {"GMX_PATH": GMX_PATH,
                   "cbpNumber": self.cbpNumber}
        
        log = "out.log"
        err = "out.err"
        self.run(SETUP_TEMPLATE, inserts, ROOT_DIR, log, err)
        self.run(GPP_TEMPLATE, inserts, join(ROOT_DIR, self.cbpNumber), log, err)
        
        self.run(MDRUN_TEMPLATE, inserts, join(ROOT_DIR, self.cbpNumber), log, err)
    
    def run(self, argString, inserts, workDir, logPath, errPath):
        
        argString = argString.format(**inserts)
        argList = argString.split()
        
        logFile = open(join(workDir, logPath),"a")
        errFile = open(join(workDir, errPath),"a")
        try:
            #print argString
            subprocess.Popen(argList, cwd=workDir, stdout=logFile, stderr=errFile).wait()
        except:
            print "Subprocess terminated with error: {0}".format(format_exc())
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

def main():
    cBPDeposition = CBPDeposition("300")
    
    print cBPDeposition.hasReachedLayer()
    
    insertHeight = cBPDeposition.getInsertHeight()
    
    xPos, yPos = cBPDeposition.getRandomPosXY()
    
    cBPModel_000 = pmx.Model(join(ROOT_DIR, CBP_PDB_000))
    cBPModel_000.translate([xPos, yPos, insertHeight])
    cBPModel_000.random_rotation()
    
    newCBP = int(cBPDeposition.cbpNumber) + 1
    cBPDeposition.model.insert_residue(newCBP, cBPModel_000.residues[0], " ")
    
    cBPDeposition.model.write(IN_PDB_FILE.format(**{"cbpNumber":newCBP}), "{0} CBP molecules".format(newCBP), 0)
    
    print cBPDeposition.hasReachedLayer()
    
    #runSystem("300")

main()