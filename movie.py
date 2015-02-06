import pymol
import pymol.experimenting
from pmx.xtc import *
from pmx import *
import os.path
import time
import subprocess

PYMOL_TEMPLATE="""\
set_view (\
     0.449250728,    0.289243400,   -0.845287681,\
    -0.884793460,    0.275103867,   -0.376111388,\
     0.123753794,    0.916872859,    0.379511178,\
    -0.000082850,   -0.000003487, -343.202331543,\
    44.792270660,   42.472671509,   34.783828735,\
   265.900756836,  420.503906250,  -20.000000000 )
    load {0}
    ray
    png {1}
"""

fn = 'init.pdb'
fn_xtc = 'md.xtc'

# pnx requires having set up the GMX_DLL variable pointing to gromacs shared libraries
# GMX_DLL="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"

# Read trajectory in
m = Model(fn)
trj = Trajectory(fn_xtc)

# Iterate over frames
tmp_base_fn = "tmp/" + fn_xtc[:-4]
for i, frame in enumerate(trj):
    trj.update( m )
    tmp_fn = tmp_base_fn + str(i) + ".pdb"
    m.write(tmp_fn)
    with open('pml/{0}.pml'.format(i), 'w') as fp:
      fp.write(PYMOL_TEMPLATE.format(tmp_fn, "png/{0:0>4d}.png".format(i) ))
    subprocess.Popen("pymol -qc -n pml/{0}.pml".format(i).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()


# Then make a movie with the pngs ...
# Source : https://trac.ffmpeg.org/wiki/Encode/H.264
# Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702
subprocess.Popen("ffmpeg -i png/%04d.png md.mp4".split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
