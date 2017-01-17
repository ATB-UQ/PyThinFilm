

# `Ir (ppy)_2 acac`

Delta enantiomer: IPAR

Lambda enantiomer: IPAS

Also `IPAR_oldcharges.itp` and `IPAS_oldcharges.itp` which have ppy charges
matching `IPR.itp` and `bIPS.itp` from and older ATB output.

The Lambda enantiomer crystal structure was obtained from:

Sergey Lamansky; Peter Djurovich; Drew Murphy; Feras Abdel-Razzaq;
Raymond Kwong; Irina Tsyba; Manfred Bortz; Becky Mui; Robert Bau; Mark E.
Thompson, Inorganic Chemistry (2001). 40, 1704-1711.

The cif file was converted to a pdb using the "mercury" program. This pdb
was used as an input to the ATB, with Ir replaced by Fe, in order to create
a template itp file. This ITP file was missing charges and charge groups.
Charges on the ligands were assigned using the charges on the equivalent
atoms of ppy (molid 6711) and a symmeterized protonated acac (molid 35990).


