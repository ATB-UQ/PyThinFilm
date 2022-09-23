.. _Validating Datasets:

Validating Datasets
===================

After your datasets have been prepared, csc-tools can be used to validate the contents of dataset directories to ensure that the datasets are suitable for contribution to the ACSC.

This can be done using the :code:`csct validate` subcommand.

To use the command to validate all properties of a dataset, execute the following command:

.. code-block:: console

    csct validate /path/to/datasets /path/to/more/datasets

.. note::
    Any number of dataset paths can be supplied to the validator.  Any supplied paths will be scanned recusrively for directories containing :code:`atbrepo.yaml` files, which will be be identified by csc-tools as simulation runs.  For convenience, ACSC recommends storing datasets in need of validation in a single parent directory.

This command will check the metadata files, directory structure, and files within subdirectories against the requirements outlined in :ref:`Structuring Datasets` and :ref:`Structuring Metadata` and inform the user of any problems.

.. note::
    csc-tools supports validation of individual dataset properties.  for a full list of options, execute the command :code:`csct validate --help`