import subprocess
import os
from io import StringIO
from os.path import join
from pathlib import Path
from traceback import format_exc
import pmx
import random
import logging
import math
import yaml
import jinja2
from time import time
import numpy as np

from PyThinFilm.helpers import res_highest_z, remove_residues, recursive_correct_paths, get_mass_dict, group_residues, \
    basename_remove_suffix, foldnorm_mean, calc_exclusions

from PyThinFilm.common import GPP_TEMPLATE, GROMPP, MDRUN, MDRUN_TEMPLATE, \
    MDRUN_TEMPLATE_GPU, K_B, ROOT_DIRS, DEFAULT_SETTING, TEMPLATE_DIR


class Deposition(object):

    def __init__(self, name, user_run_config, n_cores, debug):
        self.run_ID = 0
        self.name = name
        self.debug = debug
        self.n_cores = n_cores

        # load default run config
        with open(DEFAULT_SETTING) as fh:
            self.run_config = yaml.safe_load(fh)

        if not isinstance(user_run_config, dict):
            # load user-supplied run config file
            with open(user_run_config) as fh:
                user_run_config_dict = yaml.safe_load(fh)
        else:
            user_run_config_dict = user_run_config
        self.run_config.update(user_run_config_dict)

        self.rootdir = os.path.abspath(self.run_config["work_directory"])
        # Create the 'work_directory' if it doesn't exist
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)
            for d in ROOT_DIRS:
                os.makedirs(os.path.join(self.rootdir, d))

        # initialise logging temporarily for run_ID=0
        self.setup_logging()

        # Convert file paths in run_config to be absolute where files in resources directory are used
        recursive_correct_paths(self.run_config)

        self.last_run_ID = self.get_latest_run_ID()
        if self.last_run_ID > 0 and self.has_run_failed(self.last_run_ID):
            logging.warning("Last run {} failed, associated files will be moved out of the way".format(self.last_run_ID))
            self.delete_run(self.last_run_ID)
            self.run_ID = self.last_run_ID - 1
            if self.last_run_ID > 0 and self.has_run_failed(self.get_latest_run_ID()):
                raise Exception("The previous 2 runs were found to have not completed successfully: {} and {}".format(
                    self.last_run_ID, self.last_run_ID - 1)
                )
        else:
            self.run_ID = self.last_run_ID + 1

        self.setup_logging()

        self.gmx_executable = self.run_config['gmx_executable'] \
            if self.n_cores == 1 \
            else self.run_config['gmx_executable_mpi']

        if self.run_ID == 1:
            if "initial_config_file" in self.run_config:
                configuration_path = self.run_config["initial_config_file"]
            else:
                configuration_path = self.run_config["substrate"]["pdb_file"]
        else:
            configuration_path = self.filename("final-coordinates", "gro", prev_run=True)

        self.model = pmx.Model(str(configuration_path))
        self.model.a2nm()
        # calc_exclusions(self.model)
        # raise
        self.sampling_mixture = self.run_config["mixture"]
        self.mixture = {}

        self.new_residues = []
        self.insertion_height = None

        self.mdp_template_file = Path(self.run_config['mdp_template'])
        self.topo_template_file = Path(self.run_config['topo_template'])
        self.init_gromacs_templates()

        self.insertions_per_run = self.run_config["insertions_per_run"]

        self.init_deposition_steps()

        self.remove_top_molecule = self.run_config["remove_top_molecule"] \
            if "remove_top_molecule" in self.run_config else 0
        self.solvent_name = self.run_config["solvent_name"] \
            if "solvent_name" in self.run_config else None

    def init_gromacs_templates(self):
        if not Path(self.run_config['mdp_template']).exists():
            self.mdp_template_file = TEMPLATE_DIR / self.run_config['mdp_template']
            if not self.mdp_template_file.exists():
                raise Exception("Cannot find required file: {}. A GROMACS mdp file template is required.".format(
                    self.run_config['mdp_template']))
        if not Path(self.run_config['topo_template']).exists():
            self.topo_template_file = TEMPLATE_DIR / self.run_config['topo_template']
            if not self.topo_template_file.exists():
                raise Exception("Cannot find required file: {}. A GROMACS topo file template is required.".format(
                    self.run_config['topo_template']))

    def run_gromacs_simulation(self):

        arg_values = {"GMX_EXEC": self.gmx_executable,
                      "MDP_FILE": self.filename("control", "mdp"),
                      "tpr": self.filename("tpr", "tpr"),
                      "xtc": self.filename("trajectory", "xtc"),
                      "top": self.filename("topology", "top"),
                      "log": self.filename("log", "log"),
                      "edr": self.filename("energy", "edr"),
                      "cpo": self.filename("checkpoint", "cpt"),
                      "initial": self.filename("input-coordinates", "gro"),
                      "final": self.filename("final-coordinates", "gro"),
                      "restraints": self.filename("restraints", "gro"),
                      "run_ID": self.run_ID,
                      "mdrun": MDRUN,
                      "grompp": GROMPP,
                      }

        if "use_gpu" in self.run_config and self.run_config["use_gpu"]:
            mdrun_template = MDRUN_TEMPLATE_GPU
        else:
            mdrun_template = MDRUN_TEMPLATE

        if self.n_cores > 1:
            mdrun_template = "{} {}".format(self.run_config["mpi_template"], mdrun_template)
            arg_values["n_cores"] = self.n_cores

        # run grompp
        self.run_subprocess(GPP_TEMPLATE, arg_values)

        # run the md
        self.run_subprocess(mdrun_template, arg_values)

        configuration_path = self.filename("final-coordinates", "gro")

        if not os.path.exists(configuration_path):
            logging.error("MD run did not produce expected output file")

        # now update model after run
        self.model = pmx.Model(str(configuration_path))
        self.model.a2nm()

    def init_gromacs_simulation(self):

        if not os.path.exists(self.rootdir):
            os.mkdir(self.rootdir)

        with open(self.topo_template_file) as fh:
            top_template = jinja2.Template(fh.read())

        with open(self.mdp_template_file) as fh:
            mdp_template = jinja2.Template(fh.read())
        residue_list = list(map(lambda x:x.resname, self.model.residues))
        bulk_itps = [molecule["itp_file"] for molecule in self.mixture.values() if molecule["count"] > 0]
        substrate_itps = [self.run_config["substrate"]["itp_file"]] if "substrate" in self.run_config else []
        itp_files_include = ['#include "{}"'.format(f) for f in bulk_itps + substrate_itps]
        with open(self.filename("topology", "top"), "w") as fh:
            fh.write(
                top_template.render(
                    itp_files_include="\n".join(itp_files_include),
                    forcefield=self.run_config["forcefield_file"],
                    res_groups="\n".join(group_residues(residue_list)),
                )
            )
        n_steps = int(self.run_config["run_time"] / self.run_config["time_step"])
        if self.run_config["drift_velocity"] > 0.0:
            # Remove center of mass motion every 2x number of steps required to reach existing layer. Prior to reaching
            # the layer causes the molecule being deposited to be halted, particularly early in deposition simulations.
            mean_flight_time = self.run_config["insert_distance"] / self.run_config["drift_velocity"]
            nstcomm = 2*int(
                mean_flight_time / self.run_config["time_step"]
            )
            if nstcomm < n_steps:
                logging.debug(f"Center of mass motion will be removed every {mean_flight_time*2} ps (nstcomm={nstcomm})")
            if nstcomm == n_steps:
                # if nstcomm == n_steps (nsteps) then the deposited molecule's drift velocity is removed on the first step.
                nstcomm -= 1
                logging.debug(
                    f"Center of mass motion will be removed every {nstcomm * self.run_config['time_step'] } ps (nstcomm={nstcomm})")
            else:
                logging.warning(f"Center of mass motion will not be removed since run_time "
                                f"({self.run_config['run_time']} ps) is less than 2x average deposition time i.e. "
                                f"insert_distance/drift_velocity = {mean_flight_time} ps.")
        else:
            nstcomm = 100
        res_list = list(set(residue_list))
        mdp = mdp_template.render(res_list=" ".join(res_list),
                                  time_step=self.run_config["time_step"],
                                  n_steps=n_steps,
                                  nstcomm=nstcomm,
                                  temperature_list=" ".join([str(self.run_config["temperature"])]*len(res_list)),
                                  tau_t_list=" ".join([str(self.run_config["tau_t"])]*len(res_list)),
                                  )
        with open(self.filename("control", "mdp"), "w") as fh:
            fh.write(mdp)

    def run_subprocess(self, cli_template, inserts):

        arg_list_template = cli_template.split()
        arg_list = [arg.format(**inserts) for arg in arg_list_template]

        log_filepath = self.filename("stdout", "txt")
        err_filepath = self.filename("stderr", "txt")
        with open(log_filepath, "a") as log_file, open(err_filepath, "a") as err_file:
            logging.debug("    Running from: '{0}'".format(self.rootdir))
            logging.debug("    Running command: '{0}'".format(" ".join(arg_list)))
            logging.debug("    Process stderr ({}) and stdout ({})".format(err_filepath, log_filepath))
            proc = subprocess.Popen(arg_list, cwd=self.rootdir, stdout=log_file, stderr=err_file, env=os.environ)
            return_code = proc.wait()

        if 0 < return_code:
            msg = "Subprocess terminated with nonzero exit code: \n{0}\n\n{1}".format(" ".join(arg_list), format_exc())
            logging.error(msg)
            raise Exception(msg)

    def molecule_number(self):
        return len(self.model.residues)
    
    def setup_logging(self):
        logfile = self.filename("deposition-log", "txt")
        if self.debug:
            verbosity = logging.DEBUG
            format_log = '%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)'
        else:
            verbosity = logging.INFO
            format_log = '%(asctime)s - [%(levelname)s] - %(message)s'
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=logfile, level=verbosity, format=format_log, datefmt='%d-%m-%Y %H:%M:%S')

    def get_latest_run_ID(self):
        logdir = os.path.join(self.rootdir, "input-coordinates")
        if os.path.isdir(logdir):
            log_files = os.listdir(os.path.join(self.rootdir, "input-coordinates"))
            num = [ int(f.rsplit('_', 1)[-1].split('.')[0]) for f in log_files if not ".bak." in f  ]
            return max(num) if len(num) > 0 else 0
        else:
            return 0

    def init_deposition_steps(self):
        # Get the list of all the deposition steps
        self.last_run_ID = self.run_config["n_cycles"]
        self.init_mixture_residue_counts()

    def run_parameters(self):
        return {
            'run_ID': "{0}/{1}".format(self.run_ID, self.last_run_ID),
            'temperature': "{0} K".format(self.run_config["temperature"]),
            'run_time': "{0} ps".format(self.run_config["run_time"]),
            'drift_velocity': "{0} nm/ps".format(self.run_config["drift_velocity"]),
        }

    def write_restraints_file(self):
        """Copy current configuration and reset substrate positions"""
        restraints = self.model.copy()
        substrate = pmx.Model(str(self.run_config["substrate"]["pdb_file"]))
        substrate.a2nm()
        for i in range(len(substrate.atoms)):
            restraints.atoms[i].x = substrate.atoms[i].x
        restraints_path = join(self.rootdir, self.filename("restraints", "gro"))
        restraints.write(restraints_path, "restraints for run {0}".format(self.run_ID), 0)

    def filename(self, category, ext, prev_run=False, run_ID=None):
        run_ID = run_ID if run_ID is not None else self.run_ID
        run = run_ID - (1 if prev_run else 0)
        return join(self.rootdir, category, "{}_{}_{}.{}".format(self.name, category, run, ext))

    def get_insert_height(self):
        return self.max_layer_height() + self.run_config["insert_distance"]

    def get_random_pos_xy(self, z, xy_cutoff, z_cutoff):
        Lx, Ly = map(lambda x:[v for v in x], self.model.box[:2])
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

    def max_layer_height(self, excluded_res=None, density_fraction_cutoff=0.0):
        # now takes top of layer as first 1 nm bin with atom number density less than
        # `density_fraction_cutoff` (eg 0.1) the maximum
        binwidth = 1.0 # nm
        numbins = int(math.ceil(self.model.box[2][2]/binwidth))
        count = [ 0 ] * numbins
        substrate = self.run_config["substrate"]["res_name"]
        for a in self.model.atoms:
            if a.resname != substrate:
                bin = int(math.floor(a.x[2]/binwidth))
                count[bin] += 1
        filled = [ c > density_fraction_cutoff*max(count) for c in count ]
        binned_max_height = filled.index(False) * binwidth
        logging.debug("    Max layer height {0}".format(binned_max_height))
        return binned_max_height

    def zero_substrate_velocity(self):
        for atom in self.model.atoms:
            if atom.resname == self.run_config["substrate"]["res_name"]:
                atom.v=[0.0, 0.0, 0.0]

    def top_molecule(self, resname):
        zmax = 0
        id = -1
        for residue in self.model.residues:
            if resname == residue.resname:
                for atom in residue.atoms:
                    if atom.x[2] > zmax:
                        zmax = atom.x[2]
                        id = residue.id
        return id

    # A molecule as left layer if it is above a certain height (above the layer's mean height ??)  with a net z velocity
    def has_residue_left_layer(self, residue_ID, minimum_layer_height=-1e10, density_fraction_cutoff=0.0, layer_height=None):
        if residue_ID >=1 :
            res = self.model.residue(residue_ID)
        else :
            res = self.model.residues[residue_ID]
        escape_tolerance = self.run_config["escape_tolerance"]
        highest_atom_height = max([a.x[2] for a in res.atoms])
        # if net_z_velocity <= 0 or minimum_layer_height + escape_tolerance > highest_atom_height: return False
        if minimum_layer_height + escape_tolerance > highest_atom_height:
            return False
        max_layer_height = layer_height if not layer_height == None \
            else self.max_layer_height(res, density_fraction_cutoff=density_fraction_cutoff) #  Should it be mass weighted ? why bother??
        hasLeft = highest_atom_height > max_layer_height + escape_tolerance
        if hasLeft:
            logging.debug("    Highest Atom Height: {1}".format(res.id, highest_atom_height))
        return hasLeft

    def gen_initial_velocities_last_residue(self):
        # Note that center of mass motion removal will significantly affect
        # drift velocity when few molecules are present
        target_velocity = self.run_config["drift_velocity"]
        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B * self.run_config["temperature"] / atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            fmu = foldnorm_mean(0, sigma)
            atom.v[2] = -abs(random.gauss(0, sigma)) + fmu - target_velocity
            # logging.debug("Folded Normal Mean: {} nm/ps".format(fmu*10))
        logging.debug("Mean velocities (nm/ps): Vx={:.3f}, Vy={:.3f}, Vz={:.3f}".format(
            np.mean([a.v[0] for a in self.model.residues[-1].atoms]),
            np.mean([a.v[1] for a in self.model.residues[-1].atoms]),
            np.mean([a.v[2] for a in self.model.residues[-1].atoms]),
        ))

    def sample_mixture(self):
        if len(self.sampling_mixture) == 1:
            return self.sampling_mixture.values()[0]
        # exact composition no longer supported
        if "exact_composition" in self.run_config:
            raise Exception("parameter 'exact_composition' no longer supported")

        num_residues = len(self.model.residues)
        generator = random.Random(num_residues * self.run_ID * self.run_config["seed"])
        random_num = generator.random()
        ratioSum = sum([v["ratio"] for v in self.sampling_mixture.values()])
        cumulative = 0.0
        for v in self.sampling_mixture.values():
            cumulative += v["ratio"]
            if random_num < cumulative/ratioSum:
                return v
        raise Exception("If you are reading this something terrible has happened")

    def get_next_molecule(self):
        next_mol = self.sample_mixture()
        next_molecule = pmx.Model(str(next_mol["pdb_file"]))
        next_molecule.a2nm()

        self.mixture[next_mol["res_name"]].setdefault("count", 0)
        self.mixture[next_mol["res_name"]]["count"] += 1
        logging.info("    Inserting molecule: {0}".format(next_mol["res_name"]) )

        with open(next_mol["itp_file"]) as fh:
            cbp_itp_string = fh.read()
        mass_dict = get_mass_dict(cbp_itp_string)

        # set masses to 1.0 to avoid warning
        for atom in next_molecule.atoms:
            atom.m = mass_dict[atom.name]

        insert_height = self.get_insert_height()
        x_pos, y_pos = self.get_random_pos_xy(
            insert_height,
            self.run_config["insertion_xy_radius"],
            self.run_config["insertion_z_radius"],
            )

        next_molecule.translate([x_pos, y_pos, insert_height])
        next_molecule.random_rotation()

        # Define its insertion height as the highest Z value amongst its atoms
        self.insertion_height = max([a.x[2] for a in next_molecule.atoms])

        return next_molecule

    def remove_residue_with_id(self, residue_ID):
        residue = self.model.residue(residue_ID)
        self.model.remove_residue(residue)
        logging.debug("Removing residue: {0}".format(residue.resname))
        # DON'T Decrease the run ID by one!!
        #self.run_ID -= 1
        # Decrease the mixture counts
        self.mixture[residue.resname]["count"] -= 1

    def remove_residues(self, residues):
        # Decrease the mixture counts
        for residue in residues:
            self.mixture[residue.resname]["count"] -= 1
        remove_residues(self.model, residues)

    def remove_residue(self, residue):
        self.remove_residue_with_id(residue.id)

    def highest_z(self):
        substrate = self.run_config["substrate"]["res_name"]
        count = sum(1 for a in self.model.atoms if a.resname != substrate)
        if count > 0:
            z = max(a.x[2] for a in self.model.atoms if a.resname != substrate)
        else:
            z = 0
        logging.debug("    highest z = {0} nm".format(z))
        return z

    def resize_box(self):
        """Adjust z box dimension to maintain specified overhead void size"""
        new_Lz = math.ceil(self.highest_z()) + self.run_config["overhead_void_space"]
        if round(new_Lz) == round(self.model.box[2][2]):
            return
        old_Lz = self.model.box[2][2]
        # half_old_Lz = 0.5*old_Lz
        self.model.box[2][2] = new_Lz

        # for atom in self.model.atoms:
        #     if atom.resname == self.run_config["substrate"]["res_name"] and atom.x[2] > half_old_Lz:
        #         atom.x[2] -= old_Lz
        logging.debug("    Box size changed from z={0} nm to z={1} nm".format(old_Lz, self.model.box[2][2]))

    def write_init_configuration(self):
        self.model.write(self.filename("input-coordinates", "gro"),
                         "{0} deposited molecules".format(self.run_ID),
                         0)

    def delete_run(self, run_ID=None):
        run_ID = run_ID if run_ID is not None else self.run_ID
        for dir in ROOT_DIRS:
            if os.path.isdir(os.path.join(self.rootdir, dir)):
                files = os.listdir(os.path.join(self.rootdir, dir))
                idstr = "_{}".format(run_ID)
                for file in files:
                    if file.rsplit(".", 1)[0].endswith(idstr):
                        filepath = os.path.join(self.rootdir, dir, file)
                        os.rename(filepath, filepath+".bak."+str(time()))

    def has_run_failed(self, run_ID=None):
        log_filename = self.filename("log", "log", run_ID=run_ID)
        final_gro_filename = self.filename("final-coordinates", "gro", run_ID=run_ID)
        if not os.path.isfile(log_filename) or not os.path.isfile(final_gro_filename):
            return True
        else:
            proc = subprocess.Popen(['tail', '-n', "10", log_filename], stdout=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            lines = stdout.decode().strip().splitlines()
            try:
                last_line = lines[-1]
                return not (last_line.startswith("Finished mdrun"))
            except IndexError:
                return True

    def init_mixture_residue_counts(self):
        # Set the count to zero for all the residues in the mixtures 
        for res in self.sampling_mixture.values():
            resname = res["res_name"]
            self.mixture[resname] = res
            self.mixture[resname]["count"] = 0

        # Then, update count for residues already in the model
        for res in self.model.residues:
            # ignore the substrate residue
            if res.resname != self.run_config["substrate"]["res_name"]:
                self.mixture[res.resname]["count"] += 1
        # Log the mixture
        logging.debug('Initial mixture is: {mixture}'.format(mixture=self.mixture))

    def get_num_slabs(self):
        slab_width = self.run_config["slab_width"]
        Lz = self.model.box[2][2]
        return int(math.ceil(Lz/slab_width))

    def write_index_file(self):
        """
        write index file containing groups for thermostatting and com subtraction"
        """
        slab_width = self.run_config["slab_width"]
        num_slabs = self.get_num_slabs()
        substrate = self.run_config["substrate"]["res_name"]
        residue_z = {residue.id: res_highest_z(residue) for residue in self.model.residues} #if not residue.resname == substrate }
        slab_strings = {b:StringIO() for b in range(num_slabs)}
        slab_atoms = {b:[] for b in range(num_slabs)}
        for atom in self.model.atoms:
            if not atom.resnr in self.new_residues:
                slab_id = int(math.floor(residue_z[atom.resnr]/slab_width))
                slab_atoms[slab_id].append(atom.id)

        min_atoms_per_slab = self.run_config["min_atoms_per_slab"]
        # make sure atoms not in small groups if possible
        num_slab_atoms = sum( len(x) for x in slab_atoms.values() )
        if num_slab_atoms < min_atoms_per_slab:
            for slab_id in range(1,len(slab_atoms)):
                slab_atoms[0] += slab_atoms[slab_id]
                slab_atoms[slab_id] = []
        else:
            #keep looping until no changes being made
            consolidated = False
            while not consolidated:
                consolidated = True
                for slab_id in range(1,num_slabs):
                    atoms = slab_atoms[slab_id]
                    if len(atoms)>0 and len(atoms) < min_atoms_per_slab:
                        slab_atoms[slab_id-1] = slab_atoms[slab_id-1] + atoms
                        slab_atoms[slab_id] = []
                        consolidated = False
                
        ndx_path = self.filename("index", "ndx")
        with open(ndx_path,"w") as f:
            for slab_id in range(num_slabs):
                f.write("[ slab{} ]\n".format(slab_id))
                for atom_id in slab_atoms[slab_id]:
                    f.write("{}\n".format(atom_id))
            f.write("[ system ]\n")
            for i in range(len(self.model.atoms)):
                f.write("{}\n".format(i+1))
            f.write("[ substrate ]\n")
            for residue in self.model.residues:
                if residue.resname == substrate:
                    for atom in residue.atoms:
                        f.write("{}\n".format(atom.id))
            f.write("[ film ]\n")
            for residue in self.model.residues:
                if not residue.resname == substrate:
                    for atom in residue.atoms:
                        f.write("{}\n".format(atom.id))

    def insert_residues(self):
        # Get the next molecule and insert it into the deposition model with random position, orientation and velocity
        self.new_residues = []
        for i in range(self.insertions_per_run):
            nextMolecule = self.get_next_molecule()
            last_residue = len(self.model.residues)
            self.model.insert_residue(last_residue, nextMolecule.residues[0], " ")
            self.gen_initial_velocities_last_residue()
            self.new_residues.append(len(self.model.residues))
