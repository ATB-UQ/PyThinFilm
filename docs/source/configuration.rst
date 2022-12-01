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

GROMACS MPI Executable Name
---------------------------

:code:`gmx_executable` (default: :code:`gmx_mpi`): name of (or path to) the `GROMACS <https://www.gromacs.org/>`_ `MPI <https://www.open-mpi.org/>`_ executable which will be used by PyThinFilm.


MPI Command Template
--------------------

:code:`mpi_template` (default: "mpirun -np {n_cores}"): string specifying the form of the mpirun call.

GROMPP Command Template
-----------------------

:code:`grompp_template` (default: "{GMX_EXEC} grompp -maxwarn 2 -f {mdp} -c {initial} -r {restraints} -p {top} -o {tpr}"): string specifying the form of the grompp call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`mdp` is the input .mdp parameter file for the run, :code:`initial` is the input coordinate file, :code:`restraints` is the input restraints file, :code:`top` is the input topology file, and:code:`tpr` is the output .tpr file.

MDRUN Command Template
----------------------

:code:`mdrun_template` (default: "{GMX_EXEC} mdrun -s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final}"): string specifying the form of the mdrun call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`tpr` is the input .tpr file for the run, :code:`xtc` is the output .xtc trajectory file, :code:`edr` is the output .edr energy file, :code:`log` is the output .log file, :code:`cpo` is the output .cpo file, and :code:`final` is the output file containing the final coordinates of the run.  

Force Field File
----------------

:code:`forcefield_file` (default: gromos54a7_atb.ff/forcefield.itp): force field file to be used in the simulation.

Parameter Template
------------

:code:`mdp_template` (default: deposition.mdp.epy): Parameter file template to be used in the simulation, in .epy format.

Topology Template
-----------------

:code:`topo_template` (default: topo.top.epy): Topology file template to be used in the simulation, in .epy format. 

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

:code:`use_self` (default: True):

  insert:
    # Flag to enable/disable insertion.
    # Useful to quickly toggle without commenting everything out.
    enabled: False   # (optional - default true)

    # Set true to use own geometry (between input_min_z and input_max_z) to find new layers
    use_self: true # (optional - default false if unspecified)

    # System to source inserted layer from. Can be left out if using use_self option.
    # Will fallback to `initial_structure_file` if unset and use_self not set,
    # but this generates a warning since the mixture ratio can drift if
    # initial_structure_file is one that hasn't been equilibrated yet (combined
    # factors of periodic replication and alternating layers of high and low
    # solute concentration)
    input_gro_file: ~

    # min and max z values of the point at which to split the main system
    insert_min_z: 45 # nm - System specific.
                     #       Should generally be just above the substrate in a region
                     #       where the structure is close to that of the bulk solution.

    insert_max_z: 55 # nm - System specific.
                     #       Should generally be below the bottom of the skin density
                     #       tail when the detected skin bottom (based on
                     #       `skin_density_thresh`) is at `min_skin_height`.

    # Insertion will be performed if the bottom of the skin is below this height.
    # Generally want this to be a bit higher than insert_max_z to make sure the
    # density gradient isn't interfered with.
    min_skin_height: 70 # nm - System specific.
                        #       Should be as low as feasible to minimise system size
                        #       while allowing space between the skin density
                        #       tail and the substrate for insertion

    # min and max z values between which the inserted layer should be sourced from
    source_min_z: 45 # nm - System specific.
                     #        The point above which molecules are valid targets
                     #        to be copied into the main system as an extra
                     #        layer of solution. If `use_self` is set, this
                     #        will likely be the same as `insert_min_z`. If
                     #        using an auxiliary system from `input_gro_file`
                     #        or `initial_structure_file`, then this should be
                     #        a point above which the structure is that of the
                     #        bulk solution.

    source_max_z: 60 # nm - System specific.
                     #       The point below which molecules are valid targets
                     #       to be copied into the main system as an extra
                     #       layer of solution. See above.
                     #       Note that using a smaller `source_max_z -
                     #       source_min_z` will make selection of the inserted
                     #       layer slightly faster, but provide fewer options
                     #       to choose from.

    # Optional. Insertion will only be performed if layer_height -
    # bottom_of_skin is less than this value.
    # If unset, insertions will be performed every time the bottom of the skin
    # is detected below `min_skin_height`, and the film will continue to grow.
    max_skin_thickness: 20 # nm - System specific.

    insert_thickness: 10  # nm - Thickness of inserted layer. System specific.
                          #       This will depend on the size of the molecules
                          #       in the system. Larger molecules will require
                          #       a thicker layer to ensure that enough of them
                          #       do not cross the boundary (molecules that
                          #       cross the boundary are not inserted).
                          #       Smaller values will give faster run time,
                          #       since fewer atoms are being simulated.
    thickness_tol: 0.2    # Fractional tolerance for insert_thickness
                          # (e.g. 0.2 = accept layers within +/- 20% of insert_thickness)

    # Concentration above which `consecutive_bins` bins in a row will be used
    # to detect the bottom of the skin
    skin_density_thresh: 10 # solute atoms per nm^3 - System specific.
                            #     This will depend on the concentration
                            #     of the solute in the bulk region, and the
                            #     density of the dried film.

    consecutive_bins: 8     # Fewer consecutive bins allows detection of a
                            # thinner skin, but makes that detection less
                            # reliable. Since reliability is important for
                            # layer insertion to avoid inserting too early,
                            # more bins are generally better here so long as
                            # they don't exceed the skin thickness.

    # Maximum density of solute atoms in a slab that could be selected for
    # splitting the system. Two consecutive slabs below this density are
    # searched for, and the plane between them is where the split occurs.
    # Molecules that cross the plane are deleted, so this number can be used to
    # avoid deleting too many solute molecules.
    max_solute_density: 15 # Atoms per nm^3 - System specific.
                            #     This will depend on the concentration
                            #     of the solute in the bulk region.

    # Strategy to use when choosing a layer to insert.
    # Options are:
    #  * 'best':     Choose the layer with a height closest to insert_thickness.
    #  * 'weighted': Randomly choose a layer with a higher weighting for those
    #                 that are closer in height to insert_thickness.
    #  * 'random':   Randomly choose a layer with equal weighting.
    strategy: weighted # (optional - default 'weighted')

    # Void space to leave between existing system and inserted layer.
    # Added both above and below inserted layer.
    # Should be large enough to account for ~max. van der Waals radius.
    # Will default to 0.15 with a warning if unset
    extra_space: 0.15 # nm

    # Set true to abort mdrun and exit if insertion fails
    exit_on_failure: false

  solvent_delete:
    # As above, useful for convenient toggling
    enabled: False

    # Solute concentration above which `consecutive_bins` bins in a row will be
    # used to determine the top of the slab to delete solvent molecules from.
    # Could be different to layer insertion skin_densith_thresh to allow fine-tuning.
    # (e.g. may want a slightly larger value to remove solvent from the lower
    # portion of the region with a solute density gradient, or a much larger
    # value later in the simulation to help remove the last solvent molecules)
    density_thresh: 20  # solute atoms per nm^3 - System specific.
                        #   This will depend on the concentration of the solute
                        #   in the bulk region, the density of the dried
                        #   film, and the typical gradient of the solute density
                        #   below the skin. A value should be chosen to give a
                        #   point slightly below the top of the solute density
                        #   gradient, so that solvent molecules are removed from
                        #   below that point and do not cause potential
                        #   percolation pathways through the skin to collapse
                        #   (if they exist).
                        #
                        #   Towards the end of the simulation, for the purpose
                        #   of final drying, a large value can be used so that
                        #   the layer height is used as the skin z value. This
                        #   can be combined with a large value of `slab_height`
                        #   and `slab_lower_limit: 0` to randomly delete
                        #   solvent molecules from anywhere in the system.

    consecutive_bins: 4 # Fewer consecutive bins allows density_thresh to be
                        # detected with a thinner skin, but makes detection
                        # less reliable.

    #  Height of slab to randomly remove solvent molecules from
    slab_height: 20 # nm - System specific.
                    #       This will depend on `density_thresh`, and should be
                    #       chosen so that solvent molecules are removed from
                    #       the section of the upper section of the solute
                    #       density gradient below the skin.

    # Minimum z value below which solvent molecules should not be deleted.
    # If the bottom of the slab for deletion is below this point, it will be
    # truncated, and the number of deleted molecules will be adjusted to
    # maintain an equivalent density of deleted molecules.
    slab_lower_limit: 5 # nm - System specific.
                        #       This should initially be chosen as the point
                        #       below which the structure of the solution is
                        #       influenced by the substrate. Towards the end of
                        #       a simulation for finaly drying, it may be
                        #       chosen as 0 to enable deletion from anywhere in
                        #       the system.

    # Number of solvent molecules to delete
    number: 10 # System specific.
               #  This will depend on the size of the solvent molecules, the
               #  `slab_height`, and the x,y dimensions of the system.

    # Minimum distance between chosen molecules
    min_separation: 5 # nm - System specific.
                      #       This should be chosen to prevent nearby solvent
                      #       molecules from being deleted at the same time as
                      #       each other.

    # Set true to abort mdrun and exit if no candidates for deletion.
    # Useful to know when to begin the next stage of a simulation.
    exit_on_impossible: true