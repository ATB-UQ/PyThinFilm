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
from time import time
from threading import Timer

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

GPP_TEMPLATE = "{{GMX_PATH}}{{grompp}} -f {{MDP_FILE}} -c {struct} -r {struct} -p topo.top -o md.tpr".format(struct=IN_STRUCT_FILE)


GROMPP = "grompp_d"
TPBCONV = "tpbconv_d"
MDRUN = "mdrun_d"

MDRUN_TEMPLATE = "{{mpiRun}}{{GMX_PATH}}{{mdrun}} -pd -s md.tpr -deffnm md -c {struct}".format(struct=OUT_STRUCT_FILE)

MDRUN_TEMPLATE_GPU = "{{mpiRun}}{{GMX_PATH}}{{mdrun}} -dlb no {{domain_decomposition}} -nstlist {{neighborUpdate}} -ntomp 1 -s md.tpr -deffnm md -c {struct}".format(struct=OUT_STRUCT_FILE)


K_B = 0.00831451 #kJ / (mol K)

DEFAULT_PARAMETERS = {"lincs_order": 4,
              "lincs_iterations":1}

class Deposition(object):

    def __init__(self, runConfigFile,
            starting_deposition_number=None,
            max_cores = None): # use runConfigFile by default
        self.runConfig = yaml.load(open(runConfigFile))
        # Convert template paths in runConfig to be absolute
        recursiveCorrectPaths(self.runConfig, PROJECT_DIR)

        self.rootdir = os.path.abspath(self.runConfig["work_directory"])
        # Create the 'work_directory' if it doesn't exist
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)

        self.run_ID = self.get_latest_run_ID()
        while self.last_run_failed():
            logging.debug("    Running failed. Deleting: '{0}'".format(self.run_ID))
            self.delete_run()
            self.run_ID = self.get_latest_run_ID()

        self.first_run_ID = self.run_ID
        self.start_time = time() 

        self.disambiguate_run_config()

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

        self.sampling_mixture = self.runConfig["mixture"]
        self.mixture = {}

        self.mdp_template_file = self.runConfig['template']
        self.mdp_file = BASENAME_REMOVE_SUFFIX(self.mdp_template_file)

        self.insertions_per_run = self.runConfig["insertions_per_run"]

        self.initDepositionSteps()

        self.log = "out.log"
        self.err = "out.err"
        self.timeout = 60*60*self.runConfig["timeout"] if "timeout" in self.runConfig else 1e30

	self.remove_top_molecule = self.runConfig["remove_top_molecule"] \
		if "remove_top_molecule" in self.runConfig else 0
	self.solvent_name = self.runConfig["solvent_name"] \
		if "solvent_name" in self.runConfig else None

        self.max_cores = self.runConfig["max_cores"] if max_cores == None else max_cores

    def molecule_number(self):
        return len(self.model.residues)

    def get_latest_run_ID(self):
        depositions = [int(d) for d in os.listdir(self.rootdir) if d.isdigit()]
        depositions.sort()
        return depositions[-1] if len(depositions) > 0 else 0

    def get_rundir(self):
        suffix = str(self.runConfig['development']['branch']) if ("development" in self.runConfig) else ""
        return join(self.rootdir, str(self.run_ID) + suffix)

    def initDepositionSteps(self):
        # Get the list of all the deposition steps
        self.last_run_ID = self.runConfig["final_deposition_number"]
        self.initMixtureAndResidueCounts()

    def disambiguate_run_config(self):
        self.gmx_path = self.runConfig['gmx_path'] if self.runConfig['gmx_path'] else ""

    def runParameters(self):
        parameters_dict =  { \
                'run_ID': "{0}/{1}".format(self.run_ID, self.last_run_ID),
                'temperature':    "{0} K".format(self.runConfig["temperature"]),
                'run_time':       "{0} ps".format(self.runConfig["run_time"]),
        }
        parameters_dict['drift_velocity'] = "{0} nm/ps".format(self.runConfig["drift_velocity"])
        return parameters_dict

    def updateModel(self, configurationPath):
        self.model = pmx.Model(configurationPath)
        self.model.nm2a()

    def runSystem(self):

        if "domain_decomposition" in self.runConfig:
            domain_decomposition = "-dd "+self.runConfig["domain_decomposition"]
        else:
            domain_decomposition = ""

        if "mpi_prefix" in self.runConfig:
            max_cores = self.max_cores 
            mpi_prefix = self.runConfig["mpi_prefix"] \
                    if "mpi_prefix" in self.runConfig \
                    else ""
            mpiRun = mpi_prefix.format(max_cores)
            mdrun = self.runConfig["mdrun"] 
        else:
            mdrun = self.runConfig["mdrun"] \
                        if "mdrun" in self.runConfig else MDRUN
            mpiRun = ""
        grompp = self.runConfig["grompp"] \
                    if "grompp" in self.runConfig else GROMPP

        inserts = {"GMX_PATH":   self.gmx_path,
                   "MDP_FILE": self.mdp_file,
                   "run_ID": self.run_ID,
                   "mpiRun":   mpiRun,
                   "mdrun":    mdrun, 
                   "grompp": grompp,
                   "domain_decomposition": domain_decomposition,
                   }

        use_gpu = "use_gpu" in self.runConfig and True==self.runConfig["use_gpu"]
        if use_gpu:
            time_step = self.runConfig["time_step"]
            neighbor_list_time = self.runConfig["neighbor_list_time"]
            inserts["neighborUpdate"] = int(neighbor_list_time/time_step)
            mdrun_template = MDRUN_TEMPLATE_GPU
        else:
            mdrun_template = MDRUN_TEMPLATE

        # run grompp 
        self.runGPP(inserts)

        # run the md
        self.run(mdrun_template, inserts)

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
            fh.write(
                topTemplate.render(
                    resMixture = self.mixture,
                    substrate  = self.runConfig["substrate"],
                    forcefield = self.runConfig["forcefield"], #eg "ffG54a7.itp"
                    resnameClusters = cluster(
                        map( lambda x:x.resname, self.model.residues)
                    ),
                )
            )

        with open(join(self.rundir, self.mdp_file),"w") as fh:
            resList = [res_name for res_name, res in self.mixture.items() if res["count"] > 0]
            time_step = self.runConfig["time_step"]
            neighbor_list_time = self.runConfig["neighbor_list_time"]

            # trajectory write frequency (nsxtcout)
            trajectory_steps = 0 \
                if not "trajectory_steps" in self.runConfig \
                else self.runConfig["trajectory_steps"]

            fh.write(mdpTemplate.render(resList=resList, 
                    substrate=self.runConfig["substrate"], 
                    resLength=len(resList), 
                    timeStep=time_step, 
                    numberOfSteps=int(self.runConfig["run_time"]/time_step), 
                    temperature=self.runConfig["temperature"],
                    trajectory_steps=trajectory_steps,
                    lincs_order=self.runConfig["lincs_order"] if "lincs_order" in self.runConfig else DEFAULT_PARAMETERS["lincs_order"],
                    lincs_iterations= self.runConfig["lincs_iterations"] if "lincs_iterations" in self.runConfig else DEFAULT_PARAMETERS["lincs_iterations"],
                    neighborUpdate = int(neighbor_list_time/time_step),
                    constraints=self.runConfig["constraints"],
                    )
            )

    def run(self, argString, inserts):

        argString = argString.format(**inserts)
        argList = argString.split()

        returncode = 1
        logFilepath = join(self.rundir, self.log)
        errFilepath = join(self.rundir, self.err)
        logFile = open(logFilepath, "a")
        errFile = open(errFilepath, "a")
        try:
            logging.debug("    Running from: '{0}'".format(self.rundir))
            logging.debug("    Running command: '{0}'".format(argString))
            step_start_time = time()
            # kill if talking too long
            proc = subprocess.Popen(argList, cwd=self.rundir, stdout=logFile, stderr=errFile, env=os.environ)
            timer = Timer(self.timeout, proc.kill) #kill if taking too long
            returncode = proc.wait()
        except:
            logging.error("Subprocess terminated with error: \n{0}\n\n{1}".format(argString, format_exc()))
            raise
        finally: 
            timer.cancel()
            logFile.close()
            errFile.close()
        if 0 < returncode:
            msg = "Subprocess terminated with nonzero exit code: \n{0}\n\n{1}".format(argString, format_exc())
            logging.error(msg)
            raise Exception(msg)

    def getInsertHeight(self):
        return  self.maxLayerHeight() + self.runConfig["insert_distance"]

    def getRandomPosXY(self, z, xy_cutoff, z_cutoff):
        Lx, Ly = map(lambda x:[v*10 for v in x], self.model.box[:2])
        # take the 0th and 1st element from x and y respectively as these must be rectangular boxes
        collision = True
        x, y = 0, 0
        max_attempts = 100
        num_attempts = 0
        while collision:
            if num_attempts > max_attempts:
                msg = "Failed to insert new molecule after {} attempts. Too crowded"
                msg = msg.format(max_attempts)
                logging.error(msg)
                raise Exception(msg)
                break
            x, y = random.uniform(0.0, Lx[0]), random.uniform(0.0, Ly[1])
            collision = self.collision_check(x, y, z,
                Lx[0], Ly[1], xy_cutoff, z_cutoff
            )
            num_attempts += 1
        return x, y

    def collision_check(self, x, y, z, Lx, Ly, xycut, zcut):
        xycut2 = xycut*xycut
        zcut2 = zcut*zcut
        for a in self.model.atoms:
            dz = a.x[2]-z
            dz2 = dz*dz
            if dz2 < zcut2:
                dx = a.x[0] - x
                dy = a.x[1] - y
                dx = dx-Lx if dx >  0.5*Lx else dx
                dy = dy-Ly if dy >  0.5*Ly else dy
                dx = dx+Lx if dx < -0.5*Lx else dx
                dy = dy+Ly if dy < -0.5*Ly else dy
                dxy2 = dx*dx + dy*dy
                if dxy2 < xycut2:
                    return True
        return False

    def maxLayerHeight(self, excluded_res = None, density_fraction_cutoff = 0.0):
# now takes top of layer as first 1 nm bin with atom number density less than
# `density_fraction_cutoff` (eg 0.1) the maximum
        binwidth = 10.0 # angstroms
        numbins = int(math.ceil(self.model.box[2][2]*10/binwidth))
        count = [ 0 ] * numbins #initialize bins
        substrate = self.runConfig["substrate"]["res_name"]
        for a in self.model.atoms:
            if a.resname != substrate:
                bin = int(math.floor(a.x[2]/binwidth))
                count[bin] += 1
        maxcount = max(count)
        filled = [ c > density_fraction_cutoff*maxcount for c in count ]
        maxLayerHeight = filled.index(False) * binwidth
        logging.debug("    Max layer height {0}".format(maxLayerHeight))
        return maxLayerHeight

    def zero_substrate_velocity(self):
	for atom in self.model.atoms:
            if atom.resname == self.runConfig["substrate"]["res_name"]:
		atom.v=[0.0,0.0,0.0] 

    def hasResidueReachedLayer(self, residue_ID):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        maxLayerHeight = self.maxLayerHeight(res)
        return any([a.x[2] < maxLayerHeight + self.runConfig["contact_tolerance"] for a in res.atoms])

    def top_molecule(self, resname):
	zmax=-1e10
	id = -1
	for residue in self.model.residues:
	    if resname == residue.resname:
		for atom in residue.atoms:
		    if atom.x[2] > zmax:
		        zmax=atom.x[2]
		        id=residue.id
	return id

    # A molecule as left layer if it is above a certain height (above the layer's mean height ??)  with a net z velocity
    def hasResidueLeftLayer(self, residue_ID, minimum_layer_height = -1e10, vapor_thickness = 10, density_fraction_cutoff=0.0, layer_height=None):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        escape_tolerance = self.runConfig["escape_tolerance"]


        highest_atom_height = max([a.x[2] for a in res.atoms])
        cut_off_height = self.model.box[2][2]*10 - vapor_thickness
       # if net_z_velocity <= 0 or minimum_layer_height + escape_tolerance > highest_atom_height: return False
	if minimum_layer_height + escape_tolerance > highest_atom_height and cut_off_height > highest_atom_height: return False

        #net_z_velocity = sum( a.v[2] for a in res.atoms ) / len(res.atoms)
        maxLayerHeight = layer_height if not layer_height == None else self.maxLayerHeight(res, density_fraction_cutoff=density_fraction_cutoff) #  Should it be mass weighted ? why bother??
        hasLeft = highest_atom_height > maxLayerHeight + escape_tolerance or highest_atom_height > cut_off_height
        if hasLeft:
            #logging.debug("    Net Z velocity for residue {0}: {1}; Highest Atom Height: {2}".format(res.id, net_z_velocity, highest_atom_height))
            logging.debug("    Highest Atom Height: {1}".format(res.id, highest_atom_height))
        return hasLeft

    def genInitialVelocitiesLastResidue(self):

        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B*self.runConfig["temperature"]/atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            atom.v[2] = -abs(random.gauss(self.runConfig["drift_velocity"], sigma))

    def sampleMixture(self):
        if len(self.sampling_mixture) == 1:
            return self.sampling_mixture.values()[0]
        # exact composition no longer supported
        if "exact_composition" in self.runConfig:
            raise Exception("parameter 'exact_composition' no longer supported")

        num_residues = len(self.model.residues)
        generator = random.Random(num_residues*self.run_ID*self.runConfig["seed"])
        random_num = generator.random()
        ratioSum = sum([v["ratio"] for v in self.sampling_mixture.values()])
        cumulative = 0.0
        for v in self.sampling_mixture.values():
            cumulative += v["ratio"]
            if random_num < cumulative/ratioSum:
                return v
        raise Exception("If you are reading this something terrible has happened")

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
        xPos, yPos = self.getRandomPosXY(
            insertHeight,
            self.runConfig["insertion_xy_radius"],
            self.runConfig["insertion_z_radius"],
            )

        nextMolecule.translate([xPos, yPos, insertHeight])
        nextMolecule.random_rotation()

        # Define its insertion height as the highest Z value amongst its atoms
        self.insertionHeight = max([ a.x[2] for a in nextMolecule.atoms])

        return nextMolecule

    def removeResidueWithID(self, residue_ID):
        residue = self.model.residue(residue_ID)
        self.model.remove_residue(residue)
        logging.debug("Removing residue: {0}".format(residue.resname))
        # DON'T Decrease the run ID by one!!
        #self.run_ID -= 1
        # Decrease the mixture counts
        self.mixture[residue.resname]["count"] -= 1

    def removeResidues(self, residues):
        # Decrease the mixture counts
        for residue in residues:
            self.mixture[residue.resname]["count"] -= 1
        remove_residues(self.model, residues)

    def removeResidue(self, residue):
        self.removeResidueWithID(residue.id)

    def highest_z(self):
        substrate = self.runConfig["substrate"]["res_name"]
        count = sum(1 for a in self.model.atoms if a.resname != substrate)
        if count > 0:
            z = max(a.x[2] for a in self.model.atoms if a.resname != substrate)
        else:
            z=0
        halfheight = 0.5*self.model.box[2][2]*10
        logging.debug("    highest z = {0}".format(z))
        return z

    def resize_box(self, new_Lz):
        halfheight = 0.5*self.model.box[2][2]*10
        self.model.box[2][2] = 0.1 * (new_Lz)

        logging.debug("    box size changed from {0} to {1}".format(halfheight*2, 10*self.model.box[2][2]))

    def writeInitConfiguration(self):
        updatedPDBPath = join(self.rundir, IN_STRUCT_FILE)
        self.model.write(updatedPDBPath, "{0} deposited molecules".format(self.run_ID), 0)

    def delete_run(self):
        run_dir = os.path.join(self.rootdir ,str(self.run_ID))
        os.rename(run_dir, run_dir+".bak."+str(time()))

    def last_run_failed(self):
        log_filename = os.path.join(self.rootdir ,str(self.run_ID), "md.log")
        if not os.path.isfile(log_filename):
            return False
        else:
            proc = subprocess.Popen(['tail', '-n', "1", log_filename], stdout=subprocess.PIPE)
            lines = proc.stdout.readlines()
            last_line = lines[-1]
            return not (last_line.startswith("Finished mdrun"))

    def initMixtureAndResidueCounts(self):
        # Set the count to zero for all the residues in the mixtures 
        for res in self.sampling_mixture.values():
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

#this is needed because removing one by one is very slow in large systems
# (one by one takes 4s for 1.8 million atoms)
#modified from pmx code
def remove_residues(model,residues):
    logging.debug("num  residue: {0} num chains {1}".format(len(model.chains), len(model.residues)))
    assert len(model.chains) == 1
    chain = model.chains[0]
    for residue in residues:
        logging.debug("Removing residue: {0}".format(residue.resname))
        idx = chain.residues.index(residue)
        try:
            midx = chain.model.residues.index(residue)
        except:
            midx = -1
        del chain.residues[idx]
        del chain.model.residues[midx]
    model.atoms = []
    for r in model.residues:
        for atom in r.atoms:
            model.atoms.append(atom)
    chain.atoms = []
    for r in chain.residues:
        for atom in r.atoms:
            chain.atoms.append(atom)
    model.renumber_atoms()
    model.renumber_residues()
    chain.make_residue_tree()

def recursiveCorrectPaths(node, root_dir):
    for key, value in node.items():
        if isinstance(value, dict):
            recursiveCorrectPaths(value, root_dir)
        elif isinstance(value, list):
            for item in value:
                recursiveCorrectPaths(item, root_dir) 
        elif "template" in key:
            node[key] = join(root_dir, value)
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

def runDeposition(runConfigFile, starting_deposition_number=None,
        continuation=False, remove_bounce=False, remove_leaving_layer=True, max_cores = None,
        debug=DEBUG):
    if debug:
        verbosity = logging.DEBUG
        format_log = '%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)'
    else:
        verbosity = logging.INFO
        format_log = '%(asctime)s - [%(levelname)s] - %(message)s'
    logging.basicConfig(level=verbosity, format=format_log, datefmt='%d-%m-%Y %H:%M:%S')

    deposition = Deposition(runConfigFile,
            starting_deposition_number=starting_deposition_number,
            max_cores = max_cores,
    )
    if deposition.run_ID == deposition.last_run_ID:
        logging.error("No more depositions to run")
        raise Exception("No more depositions to run")


    walltime = 0.0
    starttime = time()
    while deposition.run_ID < deposition.last_run_ID:
        # Increment run ID
        deposition.run_ID += 1
        logging.debug("increment run_ID to '{0}'".format(deposition.run_ID))

        # Update the deposition step in case we enter a new deposition phase

        initial_layer_height = deposition.maxLayerHeight(
            density_fraction_cutoff = deposition.runConfig["density_fraction_cutoff"],
        )

        logging.debug("run_ID is '{0}' after maxLayerHeight".format(deposition.run_ID))

#NOTE THIS HAS BEEN MOVED FORWARD TO BETTER HANDLE SOLUTION DEPOSITION
	for i in range(deposition.remove_top_molecule) :
	    residue_id=deposition.top_molecule(deposition.solvent_name)		
	    logging.info("removing top molecule {0}".format(residue_id))
	    deposition.removeResidueWithID(residue_id)

        highest = deposition.highest_z()
        overhead_void_space = deposition.runConfig["overhead_void_space"]
        new_Lz = highest + overhead_void_space
        deposition.resize_box(new_Lz)

        if remove_leaving_layer :
            # Iterate over the residues and remove the ones that left the layer
            layer_height = deposition.maxLayerHeight(density_fraction_cutoff = deposition.runConfig['density_fraction_cutoff'])
	    logging.info("Layer height {0}".format(layer_height))
            leaving = [ residue for residue in deposition.model.residues[1:] \
                       if deposition.hasResidueLeftLayer(
                               residue.id,
                               minimum_layer_height = initial_layer_height,
                               layer_height=layer_height,
                               density_fraction_cutoff = deposition.runConfig["density_fraction_cutoff"],
                       )
                       ]
            deposition.removeResidues(leaving)
        logging.info("[DEPOSITION] Running deposition with parameters: {parameters_dict}".format(rid=deposition.run_ID, parameters_dict=deposition.runParameters()))
        # Get the next molecule and insert it into the deposition model with random position, orientation and velocity
        for i in range(deposition.insertions_per_run):
            nextMolecule = deposition.getNextMolecule()
            last_residue = len(deposition.model.residues)
            deposition.model.insert_residue(last_residue, nextMolecule.residues[0], " ")
            deposition.genInitialVelocitiesLastResidue()

        # create run directory and run setup make file
        deposition.runSetup()
        
        deposition.zero_substrate_velocity()

        # Write updated model to run directory
        deposition.writeInitConfiguration()

        actualMixture = ",".join([" {0}:{1}".format(r["res_name"], r["count"]) for r in deposition.mixture.values()])
        # Do first Run
        logging.info("    Current mixture is: {0}".format(actualMixture))
        deposition.runSystem()

        walltime = time() - starttime

    if deposition.run_ID >= deposition.last_run_ID:
        logging.info("Finished deposition of {0} molecules".format(deposition.last_run_ID))
    else:
        logging.info("Allowed walltime exceeded after run_id {}".format(deposition.run_ID))


def parseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--debug', dest='debug', action='store_true')
    parser.add_argument('--start', help='{int} Provide a starting deposition number different from the one in the YAML file. Used to restart a deposition. This number corresponds to the last successful deposition. Use 0 to start from scratch.')
    parser.add_argument('--remove-mol-bouncing', dest='remove_bounce',        action='store_true')
    parser.add_argument('--remove-mol-leaving',  dest='remove_leaving_layer', action='store_true')
    parser.add_argument('--max-cores', dest='max_cores', default = None, type=int,
            help='{int} Provide the maximum number of cores to use with mpi. Overrides "max_cores" in config file.')
    args = parser.parse_args()

    runDeposition(args.input,
            starting_deposition_number=int(args.start) if args.start else None,
            remove_bounce=args.remove_bounce,
            remove_leaving_layer=args.remove_leaving_layer,
            debug=args.debug,
            max_cores = args.max_cores)

if __name__=="__main__":
    parseCommandLine()
