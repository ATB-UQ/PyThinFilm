.. _Structuring Datasets:

Structuring Datasets
====================

The ACSC uses a standard directory structure to identify files associated with molecular dynamics simulations.  Each dataset submitted to the ACSC must conform to the following directory structure:

.. code-block:: text

    └──System_Name
        ├── atbrepo.yaml
        ├── control
        |   └──System_Name_control_00001.fileextension
        ├── energy
        |   └──System_Name_energy_00001.fileextension
        ├── final-coordinates
        |   └──System_Name_final-coordinates_00001.fileextension
        ├── input-coordinates
        |   └──System_Name_input-coordinates_00001.fileextension
        ├── forcefield-files
        |   └──System_Name_forcefield-files.fileextension
        ├── log
        |   └──System_Name_log_00001.fileextension
        ├── reference-coordinates
        |   └──System_Name_reference-coordinates.fileextension
        ├── topology
        |   └──System_Name_topology.fileextension
        └── trajectory
            └──System_Name_trajectory_00001.fileextension

In all cases listed above :code:`.fileextension` should be replaced with the appropriate file extension for the given file type.  In most cases, no modification of the file extension should be necessary on the part of the user.

Files containing :code:`00001` are file types which can be sequentially numbered to indicate continuation runs.  

.. warning::
    All top level directories (control, energy, etc.) must be present, even if no files of a given type are provided.  No other files or directories, hidden or otherwise, are permitted at the top level.  Ensure that no hidden files are present with :code:`ls -a`. 

System_Name
-----------

The name of the dataset directory should be a descriptive name for the simulation run.  While this is not the name that will be displayed on the ACSC website, it is the name that will be used to generate URLs pertaining to the dataset and should be used when naming all dataset files, as outlined in the template above.  

.. warning::
    Dataset directory names must be globally unique across the entire ACSC site.  Ensure that you have chosen a sufficiently descriptive name as to be unlikely to conflict with other datasets.

atbrepo.yaml
------------

This is the metadata file for the dataset, the contents of which will be explained :ref:`in a later section <Structuring Metadata>`.

control
-------

This directory should contain control/parameter/settings files for the simulation run (or segments thereof, in the case of continuation runs) (e.g., .mdp files for GROMACS, .imd files for GROMOS, .mdin files for AMBER).

.. note::
    At least one control or log file must be provided for a dataset to be included in the ACSC database. 

energy
------

This directory should contain energy trajectory files for the simulation run (or segments thereof, in the case of continuation runs).

final-coordinates
-----------------

This directory should contain output coordinate files for the simulation run (or segments thereof, in the case of continuation runs).

input-coordinates
-----------------

This directory should contain input coordinate files for the simulation run (or segments thereof, in the case of continuation runs).

.. note::
    At least one input coordinates file must be provided for a dataset to be included in the ACSC database.

forcefield-files
----------------

This directory should contain force field modifcation files for the simulation run (e.g., .ifp files for GROMACS, .itp files for GROMOS).

log
---

This directory should contain log files for the simulation run (or segments thereof, in the case of continuation runs).

.. note::
    At least one control or log file must be provided for a dataset to be included in the ACSC database. 

reference-coordinates
---------------------

Reference coordinates for the simulation run (or other coordinate files which do not meet the criteria for input or output coordinates).

topology
--------

Topology files for the simulation run.

.. note::
    At least one topology file must be provided for a dataset to be included in the ACSC database.

trajectory
----------

Coordinate trajectory files for the simulation run (or segments thereof, in the case of continuation runs).