.. _Structuring Metadata:

Structuring Metadata
====================

As mentioned in :ref:`Structuring Datasets`, the atbrepo.yaml file in each dataset directory contains the metadata necessary to describe the dataset in the ACSC database.  This metadata must be provided by the user in `YAML <https://yaml.org/>`_ format according to the following template:

.. code-block:: yaml

    title: "This will appear as the title of the simulation on the ACSC website. Should be enclosed in quotation marks."
    notes: "This will appear as a description of the simulation on the ACSC website. Should be enclosed in quotation marks.  If the data is related to a publication, the DOI of the publication can also be included in this field."
    program: The name of the program used to carry out the simulation run. Currently supported values are AMBER, GROMACS, and GROMOS.  Only one program name per dataset should be provided.  
    organization: The organization the uploading user is affiliated with, as outlined in "Contribution Prerequisites".  Currently supported values are bernhardt, chalmers, deplazes, krenske, malde, mduq, omara, smith, and yu.
    tags: Freeform text tags for the simulation, prefixed with a hyphen. Note that some tags are prefixed with "item-", as shown below
        - replicate-[number] of [total number]
        - protein-[name of protein]
        - peptide-[name of peptide]
        - lipid-[name of lipid]
        - PDB-[pdb code]
        - solvent-[name of solvent]
 
A few example metadata files are included for reference:

.. code-block:: yaml

    title: 'D6PC lipid bilayer'
    notes: 'A lipid bilayer consisting of 512 D6PC molecules. Initiated from a smaller 128-lipid DLPC equilibrated bilayer with trimmed tails. Pore spontaneously form during the simulation.'
    program: GROMACS
    organization: mduq
    tags:
        - lipid-D6PC
        - lipid bilayer
        - GROMOS
        - 54A7
        - solvent-H2O
    
.. code-block:: yaml

    title: 'Hen egg-white lysozyme'
    notes: "Hen egg-white lysozyme protein with GROMOS 54A7 in AMBER (replicate 1 of 3).  Initial structure obtained from the Protein Data Bank (PDB). PDB ID - 1AKI, URL -  https://www.rcsb.org/structure/1AKI "
    program: AMBER
    organization: mduq
    tags:
        - GROMOS
        - 54A7
        - solvent-H2O
        - replicate-1 of 3
        - PDB-1AKI
        - protein-hen egg-white lysozyme
        - protein

.. code-block:: yaml

    title: 'Alpha-helical peptide AP'
    notes: 'Alpha-helical peptide AP with GROMOS 54A7 in GROMOS (replicate 1 of 3).'
    program: GROMOS
    organization: mduq
    tags:
        - GROMOS
        - 54A7
        - solvent-H2O
        - replicate-1 of 3
        - peptide-alpha-helical peptide AP
        - peptide