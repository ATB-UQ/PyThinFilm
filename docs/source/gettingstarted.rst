Getting Started
===============

Prerequisites
-------------

Installation of PyThinFilm requires `Python <https://www.python.org//>`_ 3.6 or later, `Setuptools <https://setuptools.pypa.io/>`_ 30.3.0 or later, and a `PEP 517 <https://peps.python.org/pep-0517/>`_-compatible front end such as `pip <https://pypi.org/project/pip/>`_ 19.0 or later.

`GROMACS <https://www.gromacs.org/>`_ 2020 or later must also be installed in order to run simulations using PyThinFilm.  If `OMP <https://www.openmp.org/>`_, `MPI <https://www.open-mpi.org/>`_, or `CUDA <https://developer.nvidia.com/cuda-toolkit>`_ versions of GROMACS are desired to be used with PyThinFilm, these must also be installed.

Installation
------------

PyThinFilm can be installed using the following command:

.. code-block:: console

    pip install git+https://github.com/ATB-UQ/PyThinFilm.git

Verify that PyThinFilm has successfully installed by executing the following command: 

.. code-block:: console

    pytf --help