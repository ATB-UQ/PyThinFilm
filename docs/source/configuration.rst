.. _Configuring PyThinFilm:
Usage and Configuration
=======================

Test
====

Command-Line Usage
------------------

The :code:`pytf` command serves as the main entry point for PyThinFilm. 

Default Configuration
---------------------

Default configuration settings are specified in default_settings.yml

Any of the default configuration settings can be overridden by specifying values for them in the relevant YAML file
(http://www.yaml.org)

List of Configuration Options and Defaults
------------------------------------------

Working Directory
~~~~~~~~~~~~~~~~~

Test
*****

:code:`work_directory`: path (either absolute or relative to the PyThinFilm install directory path) where output files will be written.  default: :code:`pytf_working`

Name
~~~~

:code:`name`: a name for a particular simulation which will prefix the names of output files.  default: :code:`PyThinFilm`
