.. _Configuration Options:
List of Configuration Options
=============================

The following is a list of all configuration options which can be specified in the `YAML <https://www.yaml.org/>`_ formatted configuration file supplied as a command-line argument to run simulations, as described in :ref:`Command-Line Interface`.  Default configuration settings are specified in default_settings.yml.  Any of the default configuration settings can be overridden by specifying values for them in the configuration file supplied on the command line.

Simulation Type
---------------

:code:`simulation_type` (required): Type of simulation to perform using PyThinFilm.  Allowed values are: vacuum_deposition, solvent_evaporation, thermal_annealing

Name
----

:code:`name` (default: :code:`PyThinFilm`): a name for a particular simulation which will prefix the names of output files.  

Working Directory
-----------------

:code:`work_directory` (default: :code:`pytf_working`): path (either absolute or relative to the current working directory) where output files will be written. 

Description
-----------------------

:code:`description` (default: :code:`PyThinFilm Simulation`): description of the simulation run

GROMACS Executable Name
-----------------------

:code:`gmx_executable` (default: :code:`gmx`): name of (or path to) the `GROMACS <https://www.gromacs.org/>`_ executable which will be used by PyThinFilm.

MPI Command Template
--------------------

:code:`mpi_template` (default: ""): string specifying the form of the mpirun call e.g., "mpirun -np $NCORES".

GROMPP Command Template
-----------------------

:code:`grompp_template` (default: "{GMX_EXEC} grompp -maxwarn 2 -f {mdp} -c {initial} -r {restraints} -p {top} -o {tpr}"): string specifying the form of the grompp call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`mdp` is the input .mdp parameter file for the run, :code:`initial` is the input coordinate file, :code:`restraints` is the input restraints file, :code:`top` is the input topology file, and :code:`tpr` is the output .tpr file.

MDRUN Command Template
----------------------

:code:`mdrun_template` (default: "{GMX_EXEC} mdrun -s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final}"): string specifying the form of the mdrun call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`tpr` is the input .tpr file for the run, :code:`xtc` is the output .xtc trajectory file, :code:`edr` is the output .edr energy file, :code:`log` is the output .log file, :code:`cpo` is the output .cpo file, and :code:`final` is the output file containing the final coordinates of the run.  

Force Field File
----------------

:code:`forcefield_file` (default: gromos54a7_atb.ff/forcefield.itp): force field file to be used in the simulation.

Parameter Template
------------

:code:`mdp_template` (default: deposition.mdp.epy): parameter file template to be used in the simulation, in .epy format.

Topology Template
-----------------

:code:`topo_template` (default: topo.top.epy): topology file template to be used in the simulation, in .epy format. 

Initial Structure File
-----------------------

:code:`initial_structure_file`: starting structure file for solution evaporation and thermal annealing simulations.

Solvent Name
------------

:code:`solvent_name` (default: None): residue name of solvent molecules in simulation.

Substrate
---------

:code:`substrate` (default: None): dictionary containing the details of the substrate molecules to be used in a simulation.  Two items, :code:`res_name` and :code:`itp_file` must be specified describing the name of the substrate residues and the location of the associated ITP file, respectively.  If the :code:`initial_structure_file` configuration option is not set (in the case of vacuum deposition simulations), a third entry, :code:`pdb_file` must also be defined describing the location of the relevant pdb file.

Mixture
-------

:code:`mixture` (default: None): list of dictionaries specifying the composition of the solute mixture in solvent evaporation and vacuum deposition simulations.  For solvent evaporation and thermal annealing simulations, the dictionaries must include two entries, :code:`res_name` and :code:`itp_file` describing the residue name of the solute component and the location of the relevant itp file, respectively.  For vacuum deposition simulations, two additional entries, :code:`pdb_file` and :code:`ratio` must be supplied describing the location of the relevant pdb file and the relative abundance of the solute in the mixture to be deposited, respectively.

Number of Cycles
-----------

:code:`n_cycles` (default: 1): Number of simulation cycles to run. Ignored for thermal annealing simulations.   

Random Seed
-----------

:code:`seed` (default: 0): random number seed used in the simulation.

Temperature
-----------

:code:`temperature` (default:300): target thermostat temperature in Kelvin for the simulation.

Temperature Coupling Constant
------------------------------

:code:`tau_t` (default: 0.1): temperature coupling constant to be used in the simulation.

Temperature List
-----------------

:code:`temperature_list`: list of temperature values for thermal annealing simulations in K.  Simulations will be run for :code:`run_time` at each temperature value.

Time Step
---------

:code:`time_step` (default: 0.002): time step of the simulation in ps.

Run Time
--------

:code:`run_time` (default: 100): length of the simulation in ps.

Write Frequency
---------------

:code:`ntstout` (default: 5000): frequency (in number of steps) with which to write output files.  

Nonbonded Cutoff Distance
-------------------------

:code:`cutoff` (default: 1.4): cutoff distance for nonbonded interactions, in nm.

Dielectric Constant
---------------------

:code:`dielectric_constant` (default: 1): dielectric constant used for the reaction field in the simulation.

Insert Distance
---------------

:code:`insert_distance` (default: 5): insert distance in nm  (vacuum deposition simulations).

Escape Tolerance
----------------

:code:`escape_tolerance` (default 5.0): escape tolerance in nm (vacuum deposition and solvent evaporation simulations). 

Density Fraction Cutoff
--------------------
:code:`density_fraction_cutoff` (default: 0.0): density fraction cutoff (vacuum deposition and solvent evaporation simulations).

Overhead Void Space
----------------
:code:`overhead_void_space` (default: 10.0): overhead void space in nm (vacuum deposition and solvent evaporation simulations).

Deposition Velocity
------------------

:code:`deposition_velocity` (default: 0.0): deposition velocity in nm/ps (vacuum deposition simulations).

Insertions per Run
-------------------

:code:`insertions_per_run` (default: 1): number of insertions per run (vacuum deposition simulations).

Maximum Insertion Attempts
--------------------------

:code:`max_insertion_attempts` (default: 100): maximum number of insertion attempts (vacuum deposition simulations).


Insertion Radius (XY)
----------------------

:code:`insertion_xy_radius` (default: 2.0): insertion radius in the x-y plane in nm (vacuum deposition simulations).

Insertion Radius (Z)
--------------------

:code:`insertion_z_radius` (default: 1.0): insertion radius measured along the z axis in nm (vacuum deposition simulations).


Slab Width
----------

:code:`slab_width` (default: 100): slab width (vacuum deposition and solvent evaporation simulations).

Minimum Atoms Per Slab
----------------------

:code:`min_atoms_per_slab` (default: 1000): minimum number of atoms per slab (vacuum deposition and solvent evaporation simulations).


Number of Highest Molecules to Remove
-------------------------------------

:code:`remove_n_highest_molecules` (default: 0): number of highest molecules to remove (solvent evaporation simulations). 

Solution Acceleration Options
-----------------------------

The following options are specified under the heading :code:`solution_acceleration` and are specific to solution evaporation simulations.

Bin Size
~~~~~~~~

:code:`density_prof_bin` (default: 0.25): bin size to use when analysing density profile for skin detection and layer insertion in nm.


Insert
~~~~~~

The following options are specified under the subheading :code:`insert` and control the insertion of additional solvent layers in solvent evaporation simulations.

:code:`enabled` (default: False): controls whether additional solution layers are inserted. All other options in this category are ignored if this value is set to False.

:code:`use_self` (default: True): controls whether own geometry (between input_min_z and input_max_z) is used to find new layers.

:code:`input_gro_file` (default: ~): system to source inserted layer from.

:code:`insert_min_z` (default: 45): minimum z value of the point at which to split the main system in nm.  The point should generally be just above the substrate in a region where the structure is close to that of the bulk solution.

:code:`insert_max_z` (default: 45): maximum z value of the point at which to split the main system in nm.  Should generally be below the bottom of the skin density tail.

:code:`min_skin_height` (default: 70): insertion will be performed if the bottom of the skin is below this height in nm.

:code:`source_min_z` (default: 45): the point above which molecules are valid targets to be copied into the main system as an extra layer of solution in nm.

:code:`source_max_z` (default: 60): the point below which molecules are valid targets to be copied into the main system as an extra layer of solution in nm.

:code:`max_skin_thickness` (default: 20): insertion will only be performed if the thickness of the layer in nm is below this value.

:code:`insert_thickness:` (default: 10): thickness of inserted layer in nm. appropriate values will depend on the size of the molecules in the system, with larger molecules requiring thicker layers.  Lower values will help improve simulation speed.

:code:`thickness_tol` (default: 0.2): fractional tolerance for insert thickness.

:code:`consecutive_bins` (default: 8): number of consecutive bins which must be detected to determine the presence of a skin.  Fewer consecutive bins allows detection of a thinner skin, but makes that detection less reliable. Since reliability is important for
                            # layer insertion to avoid inserting too early,
                            # more bins are generally better here so long as
                            # they don't exceed the skin thickness.

:code:`skin_density_thresh` (default: 10): solute concentration above which :code:`consecutive_bins` bins in a row will be used to detect the bottom of the skin in units of atoms per cubic nm. This is dependent on the concentration of the solute in the bulk region and the density of the dried film.

:code:`max_solute_density` (default: 15): maximum density of solute atoms in a slab that could be selected for splitting the system in atoms per cubic nm. Two consecutive slabs below this density are searched for, and the plane between them is where the split occurs. Molecules that cross the plane are deleted, so this number can be used to avoid deleting too many solute molecules.


:code:`strategy` (default: weighted): strategy to use when choosing a solvent layer to insert.  Options are: 'best':     Choose the layer with a height closest to insert_thickness. 'weighted': Randomly choose a layer with a higher weighting for those that are closer in height to insert_thickness.  'random':   Randomly choose a layer with equal weighting.

:code:`extra space` (default: 0.15): void space to leave between existing and inserted layer in nm. This value should be large enough to account for ~max. van der Waals radius.

:code:`exit_on_failure` (default: False): set this value to true to abort mdrun and exit if insertion fails.

Solvent Deletion
~~~~~~

The following options are specified under the subheading :code:`solvent_delete` and control the deletion of non-evaporated solvent in solvent evaporation simulations.

:code:`enabled` (default: False): toggles non-evaporated solvent deletion on or off.  If set to False, all other solvent deletion configuration settings are ignored.  

:code:`density_thresh` (default: 20): solute concentration in atoms per cubic nm above which :code:`consecutive_bins` bins in a row will be used to  (e.g. may want a slightly larger value to remove solvent from the lower portion of the region with a solute density gradient, or a much larger value later in the simulation to help remove the last solvent molecules)

:code:`consecutive_bins` (default: 4): fewer consecutive bins allows :code:`density_thresh`` to be detected with a thinner skin, but makes detection less reliable.

:code:`slab_height` (default: 20) height of slab in nm. This will depend on :code:`density_thresh`, and should be chosen so that solvent molecules are removed from the section of the upper section of the solute density gradient below the skin.

:code:`slab_lower_limit` (default: 5): minimum z value in nm below which solvent molecules should not be deleted. If the bottom of the slab for deletion is below this point, it will be truncated, and the number of deleted molecules will be adjusted to maintain an equivalent density of deleted molecules.  

:code:`number` (default: 10): number of solvent molecules to delete.  Appropriate values are related to the size of the solvent molecule, :code:`slab_height`, and the x-y dimensions of the system.  

:code:`min_separation` (default: 5): minimum distance between solvent molecules selected for deletion.  This is intended to prevent nearby solvent molecules from being deleted at the same time as each other.

:code:`exit_on_impossible` (default: True): if true, aborts mdrum and exits if no candidates for deletion are available.