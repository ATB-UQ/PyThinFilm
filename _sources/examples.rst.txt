.. _Examples & Helpful Tips:
Examples & Helpful Tips
=======================

Example Configuration Files
---------------------------

PyThinFilm includes representative examples of configuration files for each of the three simulation types the package is capable of performing, along with the necessary input files to perform these simulations, given that the criteria for default GROMACS and mpirun paths outlined in :ref:`Getting Started` are met.  

These examples are intended to serve as a starting point for users wishing to set up their own simulations by modifying the values of these configuration options.  These can be found in the :code:`src/PyThinFilm/examples` directory.  Note that these examples rely on a mixture of fallback behaviour to the default configuration options described in :ref:`Configuration Options`, overriding of these default values, and specification of configuration option values not encompassed by the defaults to function correctly.  Also note that larger values of :code:`run_time` and :code:`nstout` than are used in these examples will be required to achieve good performance in production simulations.  

Managing Disk Space
-------------------

PyThinFilm can generate a large number of files if a sufficient number of rounds of vacuum deposition or solvent evaporation are performed.  This can be mitigated by setting the GROMACS :code:`GMX_MAXBACKUP` environment variable to a value of -1.  Additionally, the accumulation of unnecessary files can be managed using a script located in :code:`src/PyThinFilm/scripts/cleanup_files.py` which removes old output files after a large number of rounds of vacuum deposition or solution evaporation have been run.