.. _Getting Started:
Getting Started
===============

Prerequisites
-------------

Installation of PyThinFilm requires `Python <https://www.python.org/>`_ 3.6 or later, `Setuptools <https://setuptools.pypa.io/>`_ 30.3.0 or later, and a `PEP 517 <https://peps.python.org/pep-0517/>`_-compatible front end such as `pip <https://pypi.org/project/pip/>`_ 19.0 or later.

`GROMACS <https://www.gromacs.org/>`_ 2020 or later must also be installed in order to run simulations using PyThinFilm.  If `OMP <https://www.openmp.org/>`_, `MPI <https://www.open-mpi.org/>`_, or `CUDA <https://developer.nvidia.com/cuda-toolkit>`_ versions of GROMACS are desired to be used with PyThinFilm, these must also be installed.

Installation
------------

Using pip
~~~~~~~~~

PyThinFilm can be installed using pip by executing the following command:

.. code-block:: console

    pip install git+https://github.com/ATB-UQ/PyThinFilm.git

From Source
~~~~~~~~~~

Alternatively, PyThinFilm can be installed from source using the following commands:

.. code-block:: console

    git clone https://github.com/ATB-UQ/PyThinFilm.git
    python setup.py install

For Development
~~~~~~~~~~~~~~

Users wishing to modify or extend PyThinFilm can install the package in editable mode by first cloning the repository and entering into the PyThinFilm directory: 

.. code-block:: console

    git clone https://github.com/ATB-UQ/PyThinFilm.git
    cd PyThinFilm

And then installing the package in editable mode using either pip:

.. code-block:: console

    pip install -e .

Or directly using Python:

.. code-block:: console

    python setup.py develop

Testing PyThinFilm
~~~~~~~~~~~~~~~~~~

After selecting an installation method, verify that PyThinFilm has successfully installed by executing the following command: 

.. code-block:: console

    pytf --help