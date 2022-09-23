.. _Exporting Datasets:

Exporting Datasets
==================

After datasets have been successfully validated, csc-tools can be used to archive datasets in preparation for upload to the ACSC server.

This can be done by passing a flag to the :code:`csct validate` subcommand:

.. code-block:: console

    csct validate --export /path/to/datasets /path/to/more/datasets

.. warning::
    Datasets must have passed validation as described in :ref:`Validating Datasets` in order to be exported and the export_path configuration option must be set as described in :ref:`Configuring csc-tools`

The progress of the dataset export can be monitored in the terminal output.  After export has completed, each successfully exported dataset will result in the creation of a corresponding :code:`.tar.gz` file in :code:`export_path` containing all files and subdirectories in the dataset directory, along with the metadata file.  

.. note::
    The export process does not delete or modify the original files in the dataset paths passed to :code:`csct validate`