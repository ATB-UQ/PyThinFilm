Installing csc-tools
====================

The ACSC provides a `Python <https://www.python.org/>`_ command-line tool to aid contributors in preparing data for upload.  This tool provides a means of validating that user-supplied data is in the correct format for upload to the ACSC servers.

Requirements
------------

csc-tools requires `Python <https://www.python.org/>`_ version 3.6 or higher.  To check if your environment is configured with an appropriate version of Python, execute the following command in a terminal:


.. code-block::

    python --version

Installation of csc-tools also requires `pip <https://pypi.org/project/pip/>`_ version 10.0 or higher.  To check if your environment is configured with an appropriate version of pip, execute the following command in a terminal:

.. code-block::

    pip --version

Installation
------------

To install csc-tools, execute the following command in a terminal:

.. code-block:: console

    pip install --upgrade --force-reinstall git+https://github.com/ATB-UQ/csc-tools@master

.. note::
    The :code:`--upgrade` and :code:`--force-reinstall` flags ensure that csc-tools is reinstalled in cases where a preexisting installation of csc-tools is present.

To confirm that the csc-tools installation is working, try running the program by executing the following command in a terminal:

.. code-block:: console

    csct --version

If the version number of :code:`csct` is displayed, the installation has completed successfully and csc-tools is ready to be used.
