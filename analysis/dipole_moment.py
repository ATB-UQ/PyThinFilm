#!/usr/bin/env python

import sys
from math import sqrt

itpfile, pdbfile = sys.argv[1:]

with open(itpfile) as f: itp = f.readlines()
with open(pdbfile) as f: pdb = f.readlines()

itp = [line for line in itp if not line[0] == ";"]

charges = []

atomsection = False

for line in itp:
    if not atomsection and line == "[ atoms ]\n":
        atomsection = True
    elif atomsection:
        words = line.split()
        if len(words) >= 8:
            charges.append( float(words[6]) )
        else: #end of section
            break

pdb = [ line for line in pdb if len(line)>4 and line[:4] in ("ATOM", "HETA") ]

x = []
y = []
z = []

for line in pdb:
    words = line.split()
    x.append(float(words[5]))
    y.append(float(words[6]))
    z.append(float(words[7]))

assert len(x) == len(charges)

dx = sum(i * qi for i,qi in zip(x, charges))
dy = sum(i * qi for i,qi in zip(y, charges))
dz = sum(i * qi for i,qi in zip(z, charges))

eA = 1/0.20819434 #Debye

print("dipole: {} Debye".format( sqrt(dx**2 + dy**2 + dz**2) * eA))
print("total charge: {} e".format( sum(charges) ))



