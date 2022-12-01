List of Configuration Options
=============================

The following is a list of all configuration options which can be specified in the `YAML <https://www.yaml.org/>`_ formatted configuration file supplied as a command-line argument to run simulations, as described in :ref:`Command-Line Interface`.  Default configuration settings are specified in default_settings.yml.  Any of the default configuration settings can be overridden by specifying values for them in the configuration file supplied on the command line.

General Configuration Options
-------

The following configuration settings apply to all simulation types supported by PyThinFilm.

Simulation Type
~~~~~~~~~~~~~~~

:code:`simulation_type` (required): Type of simulation to perform using PyThinFilm.  Allowed values are: :code:`vacuum_deposition`, :code:`solution_deposition`, :code:`thermal_annealing`

Name
~~~~

:code:`name` (default: :code:`PyThinFilm`): a name for a particular simulation which will prefix the names of output files.  

Working Directory
~~~~~~~~~~~~~~~~~

:code:`work_directory` (default: :code:`pytf_working`): path (either absolute or relative to the current working directory) where output files will be written.  default: 

Description
~~~~~~~~~~~ 

:code:`description` (default: :code:`PyThinFilm Simulation`):

GROMACS Executable Name
~~~~~~~~~~~~~~~~~~~~~~~

:code:`gmx_executable` (default: :code:`gmx`): name of (or path to) the `GROMACS <https://www.gromacs.org/>`_ executable which will be used by PyThinFilm.

GROMACS MPI Executable Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:code:`gmx_executable` (default: :code:`gmx_mpi`): name of (or path to) the `GROMACS <https://www.gromacs.org/>`_ `MPI <https://www.open-mpi.org/>`_ executable which will be used by PyThinFilm.


MPI Command Template
~~~~~~~~~~~~~~~~~~~~

:code:`mpi_template` (default: "mpirun -np {n_cores}"): string specifying the form of the mpirun call.

GROMPP Command Template
~~~~~~~~~~~~~~~~~~~~~~~

:code:`grompp_template` (default: "{GMX_EXEC} grompp -maxwarn 2 -f {mdp} -c {initial} -r {restraints} -p {top} -o {tpr}"): string specifying the form of the grompp call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`mdp` is the input .mdp file for the run, :code:`initial` is the input coordinate file, :code:`restraints` is the input restraints file, :code:`top` is the input topology file, and:code:`tpr` is the output .tpr file.

MDRUN Command Template
~~~~~~~~~~~~~~~~~~~~~~

:code:`mdrun_template` (default: "{GMX_EXEC} mdrun -s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final}"): string specifying the form of the mdrun call on the command line, where :code:`GMX_EXEC` is the GROMACS executable, :code:`tpr` is the input .tpr file for the run, :code:`xtc` is the output .xtc trajectory file, :code:`edr` is the output .edr energy file, :code:`log` is the output .log file, :code:`cpo` is the output .cpo file, and :code:`final` is the output file containing the final coordinates of the run.  

Solvent Name
~~~~~~~~~~~~

:code:`solvent_name` (default: None): residue name of solvent molecules in simulation.

Substrate
~~~~~~~~~

:code:`substrate` (default: None): details of the substrate molecules to be used in a simulation.  Two fields, :code:`res_name` and

Mixture
~~~~~~~

:code:`mixture` (default: None):


Vacuum Deposition Simulations
-----------------

In addition to the general configuration options, the following additional configuration options apply to vacuum deposition simulations.

Solution Deposition Simulations
------------------------------

In addition to the general configuration options, the following additional configuration options apply to solution deposition simulations.
