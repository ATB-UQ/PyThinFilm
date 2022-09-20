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

from PyThinFilm.helpers import res_highest_z, remove_residues_faster, recursive_correct_paths, get_mass_dict, \
    group_residues, \
    foldnorm_mean, calc_exclusions

from PyThinFilm.common import GPP_TEMPLATE, GROMPP, MDRUN, MDRUN_TEMPLATE, \
    MDRUN_TEMPLATE_GPU, K_B, ROOT_DIRS, DEFAULT_SETTING, TEMPLATE_DIR, SIMULATION_TYPES, SOLVENT_EVAPORATION


class Deposition(object):

    def __init__(self, user_run_config, n_cores, debug):
        self.run_ID = 0
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
        self.name = self.run_config["name"]
        assert self.run_config["simulation_type"] in SIMULATION_TYPES, \
            f"Unexpected simulation_type ({self.run_config['simulation_type']}), " \
            f"the recognised values are: {SIMULATION_TYPES}"
        self.type = self.run_config["simulation_type"]

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
        self.density_fraction_cutoff = self.run_config["density_fraction_cutoff"]

        self.mdp_template_file = Path(self.run_config['mdp_template'])
        self.topo_template_file = Path(self.run_config['topo_template'])
        self.init_gromacs_templates()

        self.insertions_per_run = self.run_config["insertions_per_run"]
        self.escape_tolerance = self.run_config["escape_tolerance"]

        self.remove_n_highest_molecules = self.run_config["remove_n_highest_molecules"]

        self.solvent_name = self.run_config["solvent_name"]

        self.init_deposition_steps()

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
        if self.run_config["deposition_velocity"] > 0.0:
            # Remove center of mass motion every 2x number of steps required to reach existing layer. Prior to reaching
            # the layer causes the molecule being deposited to be halted, particularly early in deposition simulations.
            mean_flight_time = self.run_config["insert_distance"] / self.run_config["deposition_velocity"]
            nstcomm = 2*int(
                mean_flight_time / self.run_config["time_step"]
            )
            if nstcomm < n_steps:
                logging.debug(f"Center of mass motion will be removed every {mean_flight_time*2} ps (nstcomm={nstcomm})")
            elif nstcomm == n_steps:
                # if nstcomm == n_steps (nsteps) then the deposited molecule's deposition velocity is removed on the first step.
                nstcomm -= 1
                logging.debug(
                    f"Center of mass motion will be removed every {nstcomm * self.run_config['time_step'] } ps (nstcomm={nstcomm})")
            else:
                logging.warning(f"Center of mass motion will not be removed since run_time "
                                f"({self.run_config['run_time']} ps) is less than 2x average deposition time i.e. "
                                f"insert_distance/deposition_velocity = {mean_flight_time} ps.")
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
            num = [int(f.rsplit('_', 1)[-1].split('.')[0]) for f in log_files if not ".bak." in f]
            return max(num) if len(num) > 0 else 0
        else:
            return 0

    def init_deposition_steps(self):
        # Get the list of all the deposition steps
        self.last_run_ID = self.run_config["n_cycles"]
        self.init_mixture_residue_counts()

    def run_config_summary(self):
        return yaml.dump({
            'simulation_type': self.type,
            'run_ID': "{0}/{1}".format(self.run_ID, self.last_run_ID),
            'temperature': "{0} K".format(self.run_config["temperature"]),
            'run_time': "{0} ps".format(self.run_config["run_time"]),
            'deposition_velocity': "{0} nm/ps".format(self.run_config["deposition_velocity"]),
        })

    def write_restraints_file(self):
        """Copy current configuration and reset substrate positions"""
        logging.debug("Creating restraints file")
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

    def get_insert_height(self, layer_height=None):
        layer_height = layer_height if layer_height is not None else self.layer_height()
        return layer_height + self.run_config["insert_distance"]

    def get_random_pos_xy(self, z, xy_cutoff, z_cutoff):
        Lx = self.model.box[0][0]
        Ly = self.model.box[1][1]
        # take the 0th and 1st element from x and y respectively as these must be rectangular boxes
        collision = True
        x, y = 0, 0
        num_attempts = 0
        while collision:
            if num_attempts > self.run_config["max_insertion_attempts"]:
                msg = f"Failed to insert new molecule after {self.run_config['max_insertion_attempts']} attempts. Too crowded"
                logging.error(msg)
                raise Exception(msg)
            x, y = random.uniform(0.0, Lx), random.uniform(0.0, Ly)
            collision = self.collision_check(x, y, z, Lx, Ly, xy_cutoff, z_cutoff)
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

    def layer_height(self, excluded_resnames=None, density_fraction_cutoff=None, bin_width=0.2):
        """Find the top of the deposited layer. In some cases it is useful to define the top of the layer as the
         as first z-bin with an atom number density less than a specified fraction of the maximum density e.g.
        `density_fraction_cutoff` 0.1 corresponds to a density of 10% of the maximum (including substrate).
        """
        density_fraction_cutoff = self.density_fraction_cutoff \
            if density_fraction_cutoff is None \
            else density_fraction_cutoff
        excluded_resnames = excluded_resnames if excluded_resnames is not None else []
        binned_max_height = 0.0
        atom_counts = {}
        for residue in self.model.residues:
            if residue.resname in excluded_resnames:
                continue
            for a in residue.atoms:
                z_bin = int(math.ceil(a.x[2]/bin_width))
                atom_counts.setdefault(z_bin, 0)
                atom_counts[z_bin] += 1
        max_atom_density = max(atom_counts.values())
        bins_with_atoms = {z_bin: atom_counts
                           for z_bin, atom_counts in sorted(atom_counts.items())
                           if atom_counts > density_fraction_cutoff * max_atom_density
                           }
        if len(bins_with_atoms) == 0:
            logging.debug("    No atoms found")
        else:
            if len(bins_with_atoms) == atom_counts:
                logging.warning("    No void above deposited layer found")
            z_bin, atom_count = sorted(bins_with_atoms.items())[-1]
            binned_max_height = z_bin * bin_width
            logging.debug("    Max layer height {:.3f} nm ({} atoms/bin)".format(binned_max_height, atom_count))
        return binned_max_height

    def zero_substrate_velocity(self):
        for atom in self.model.atoms:
            if atom.resname == self.run_config["substrate"]["res_name"]:
                atom.v = [0.0, 0.0, 0.0]

    def has_residue_left_layer_solvent_evaporation(self, residue, solvent_excluded_layer_height):
        min_atom_height = np.min([a.x[2] for a in residue.atoms])
        return min_atom_height > solvent_excluded_layer_height

    def has_residue_left_layer(self, residue, layer_height):
        mean_atom_height = np.mean([a.x[2] for a in residue.atoms])
        if mean_atom_height > layer_height + self.escape_tolerance:
            logging.debug(f"    Molecule ({residue.id}) found in the gas phase, mean atom height ({mean_atom_height:.3f}) >"
                          f" layer_height ({layer_height} nm) + escape_tolerance ({self.escape_tolerance} nm)")
            return True
        else:
            return False

    def gen_initial_velocities_last_residue(self):
        # Note that center of mass motion removal will significantly affect
        # deposition velocity when few molecules are present
        target_velocity = self.run_config["deposition_velocity"]
        for atom in self.model.residues[-1].atoms:
            sigma = math.sqrt(K_B * self.run_config["temperature"] / atom.m)
            atom.v[0] = random.gauss(0.0, sigma)
            atom.v[1] = random.gauss(0.0, sigma)
            fmu = foldnorm_mean(0, sigma)
            atom.v[2] = -abs(random.gauss(0, sigma)) + fmu - target_velocity
            # logging.debug("Folded Normal Mean: {} nm/ps".format(fmu*10))

        # If the net z-velocity is +ve, which can regularly occur for low deposition velocities,
        # invert the sign of atom velocities until the mean is -ve.
        mean_z = np.mean([a.v[2] for a in self.model.residues[-1].atoms])
        if mean_z > 0:
            logging.debug(f"Mean z velocity +ve: {mean_z}")
            for a in self.model.residues[-1].atoms:
                a.v[2] = -1*a.v[2]
                mean_z = np.mean([a.v[2] for a in self.model.residues[-1].atoms])
                if mean_z < 0:
                    logging.debug(f"Mean z velocity is now -ve: {mean_z}")
                    break

        logging.debug("Mean velocities (nm/ps): Vx={:.3f}, Vy={:.3f}, Vz={:.3f}".format(
            np.mean([a.v[0] for a in self.model.residues[-1].atoms]),
            np.mean([a.v[1] for a in self.model.residues[-1].atoms]),
            np.mean([a.v[2] for a in self.model.residues[-1].atoms]),
        ))

    def sample_mixture(self):
        if len(self.sampling_mixture) == 1:
            return list(self.sampling_mixture.values())[0]
        # exact composition no longer supported
        if "exact_composition" in self.run_config:
            raise Exception("parameter 'exact_composition' no longer supported")

        num_residues = len(self.model.residues)
        generator = random.Random(num_residues * self.run_ID * self.run_config["seed"])
        random_num = generator.random()
        ratio_sum = sum([v["ratio"] for v in self.sampling_mixture.values()])
        cumulative = 0.0
        for v in self.sampling_mixture.values():
            cumulative += v["ratio"]
            if random_num < cumulative/ratio_sum:
                return v
        raise Exception("If you are reading this something terrible has happened")

    def get_next_molecule(self, layer_height=None):
        next_molecule_resname = self.sample_mixture()
        self.mixture[next_molecule_resname["res_name"]].setdefault("count", 0)
        self.mixture[next_molecule_resname["res_name"]]["count"] += 1

        next_molecule = pmx.Model(str(next_molecule_resname["pdb_file"]))
        next_molecule.a2nm()

        self.set_masses(next_molecule, next_molecule_resname)
        next_molecule.com()

        insert_height = self.get_insert_height(layer_height=layer_height)
        x_pos, y_pos = self.get_random_pos_xy(
            insert_height,
            self.run_config["insertion_xy_radius"],
            self.run_config["insertion_z_radius"],
            )
        # logging.debug(f"Initial COG of new molecule: {np.mean([a.x for a in next_molecule.atoms], axis=0)}")
        next_molecule.random_rotation()
        # logging.debug(f"Post rotation COG of new molecule: {np.mean([a.x for a in next_molecule.atoms], axis=0)}")
        next_molecule.translate([x_pos, y_pos, insert_height])
        # logging.debug(f"Post translation COG of new molecule: {np.mean([a.x for a in next_molecule.atoms], axis=0)}")
        logging.debug(f"    Inserting molecule {next_molecule_resname['res_name']}: "
                      f"{np.mean([a.x for a in next_molecule.atoms], axis=0)}")

        return next_molecule

    def set_masses(self, next_molecule, next_molecule_resname):
        with open(next_molecule_resname["itp_file"]) as fh:
            itp_string = fh.read()
        mass_dict = get_mass_dict(itp_string)
        for atom in next_molecule.atoms:
            atom.m = mass_dict[atom.name]

    def remove_gas_phase_residues(self):
        # Iterate over the residues and remove the ones that have left the layer
        gas_phase_residues = []
        if self.run_config["simulation_type"] == SOLVENT_EVAPORATION:
            excluded_resname = self.run_config["solvent_name"]
            solvent_excluded_layer_height = self.layer_height(excluded_resnames=[excluded_resname])
            logging.info(f"Layer height excluding {excluded_resname} is {solvent_excluded_layer_height} nm")
            for residue in self.model.residues:
                if residue.resname == self.run_config["solvent_name"]:
                    if self.has_residue_left_layer_solvent_evaporation(residue, solvent_excluded_layer_height):
                        gas_phase_residues.append(residue)
        else:
            # For non-solvent evaporation simulations we may also need to remove molecules from gas phase,
            # here we use a distance threshold above the existing layer (escape tolerance) to indicate whether a
            # molecule has entered the gas phase.

            # First check if molecule(s) could require removal.
            if self.highest_z() < self.escape_tolerance:
                return

            # We need to calculation layer height with a non-zero density_fraction_cutoff so that molecules that are
            # in the gas phase can be excluded from layer_height calculation.
            density_fraction_cutoff = 0.01 if self.density_fraction_cutoff == 0.0 else self.density_fraction_cutoff
            # Remove substrate so that only molecules in the layer are included in the maximum density.
            # substrate_resname = None \
            #     if self.run_config["substrate"] is None \
            #     else self.run_config["substrate"]["res_name"]
            layer_height = self.layer_height(#excluded_resnames=[substrate_resname],
                                             density_fraction_cutoff=density_fraction_cutoff)
            logging.info(f"Layer height with density threshold {100*density_fraction_cutoff}% of maximum "
                         f"is {layer_height} nm")
            for residue in self.model.residues:
                if self.has_residue_left_layer(residue, layer_height):
                    gas_phase_residues.append(residue)

        # A solvent evaporation speedup can be applied whereby a given number of molecules are
        # removed in each cycle.
        if self.remove_n_highest_molecules > 0:
            n_highest_solvent_molecules = self.n_highest_solvent_molecules()
        else:
            n_highest_solvent_molecules = []
        residues_to_remove = gas_phase_residues + n_highest_solvent_molecules
        if len(residues_to_remove) > 0:
            self.remove_residues(residues_to_remove)

    def top_molecule_n_molecules(self, resname):
        # get all z_coords and residue IDs of atoms with resname
        selected_residues = {}
        for residue in self.model.residues:
            if resname == residue.resname:
                for atom in residue.atoms:
                    selected_residues.setdefault(residue.id, []).append(atom.x[2])
        # all we need is the maximum z coord and res_id
        residue_maxz = [(max(zs), res_id) for res_id, zs in selected_residues.items()]
        # return the highest n (self.remove_n_highest_molecules) residues
        return [self.model.residues(mol[1])
                for mol in sorted(residue_maxz, reverse=True)[:self.remove_n_highest_molecules]
                ]

    def n_highest_solvent_molecules(self):
        highest_residues = self.top_molecule_n_molecules(self.solvent_name)
        logging.debug(f"Highest {self.remove_n_highest_molecules} solvent ({self.solvent_name}) "
                      f"molecules will be removed")
        return highest_residues

    def remove_residue_with_id(self, residue_id):
        residue = self.model.residue(residue_id)
        self.model.remove_residue(residue)
        logging.debug("Removing residue with ID: {0}".format(residue.id))
        # Decrease the mixture counts
        self.mixture[residue.resname]["count"] -= 1

    def remove_residues(self, residues):
        # Decrease the mixture counts
        for residue in residues:
            self.mixture[residue.resname]["count"] -= 1
        remove_residues_faster(self.model, residues)

    def remove_residue(self, residue):
        self.remove_residue_with_id(residue.id)

    def highest_z(self):
        return max(a.x[2] for a in self.model.atoms)

    def resize_box(self):
        """Adjust z box dimension to maintain specified overhead void size, minimum increment 1nm"""
        max_z = self.highest_z()
        logging.debug(f"Maximum atom height: {max_z:.3f} nm")
        new_Lz = math.ceil(max_z) + self.run_config["overhead_void_space"]
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
        """Insert molecules into the model with random position, orientation and velocity"""
        self.new_residues = []
        # Calculate layer height prior to adding any new molecules
        layer_height = self.layer_height()
        for _ in range(self.insertions_per_run):
            next_molecule = self.get_next_molecule(layer_height=layer_height)
            last_residue = len(self.model.residues)
            self.model.insert_residue(last_residue, next_molecule.residues[0], " ")
            self.gen_initial_velocities_last_residue()
            self.new_residues.append(len(self.model.residues))
        logging.info("    Current composition: " + ",".join([" {0}:{1}".format(r["res_name"], r["count"])
                                                             for r in self.mixture.values()])
                     )
