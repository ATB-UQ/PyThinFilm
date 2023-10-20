.. _Command-Line Interface:
Command-Line Interface
======================

Entry Point
-----------

The :code:`pytf` command serves as the main entry point for PyThinFilm.  

Arguments
---------

:code:`pytf` accepts a sole argument, :code:`CONFIG`, specifiying the relative or absolute path to a configuration file `YAML <http://www.yaml.org>`_, the structure of which is described in :ref:`Configuration Options`.   Example configuration files are located in :code:`src/PyThinFilm/examples` and described in :ref:`Examples & Helpful Tips`.  

.. note::
    :code:`CONFIG` must be supplied in order for :code:`pytf` to run successfully.

Command-Line Options
--------------------

:code:`pytf` can optionally be run using the :code:`--debug` flag to print debugging information, or the :code:`--help` flag to print the help menu.