.. _Configuring csc-tools:

Configuring csc-tools
=====================

Now that csc-tools has been successfully installed, configuration is required before the program can be used to validate datasets.  

Initializing the  Configuration File
------------------------------------

Before csc-tools configuration options can be set, the csc-tools configuration file must be initialized.  To do this, execute the following command:

.. code-block:: console

    csct config --init

.. warning::
    If a csc-tools configuration file already exists, this command will delete all configuration options.

Configuring the API Token
-------------------------

csc-tools must be configured with the API token generated in :ref:`Creating an ACSC API Token` in order to correctly associate contributed datasets with users.  To do this, execute the following command:

.. code-block:: console

    csct config authorization xxxxxxxxxxxxxxxxxxx

Where :code:`xxxxxxxxxxxxxxxxxxx` is the API token generated in :ref:`Creating an ACSC API Token`.  To check if the API token was successfully set, execute the following command:

.. code-block:: console

    csct config authorization

If the API token was configured correctly, this command should echo your API token.

Configuring the Export Path
---------------------------

csc-tools exports successfully validated datasets as compressed archives to a specified export path.  This export path must be configured before datasets can be prepared for upload using csc-tools.  To specify the export path for csc-tools, execute the following command:

.. code-block:: console

    csct config export_path /path/to/where/files/should/be/exported

.. note::
    Ensure that the specified export path is on a filesystem with sufficient space to hold the processed dataset archives.

To check if the export path was successfully set, execute the following command:

.. code-block:: console

    csct config export_path

If the export path was configured correctly, this command should echo your export path.