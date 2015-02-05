import pymol
import pymol.experimenting
from pmx.xtc import *

fn = 'test.pdb'

# GMX_DLL="/Users/BertrandCaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"

# Initialize a pymol session
pymol.pymol_argv = ['pymol','-qc']
pymol.finish_launching()
cmd = pymol.cmd

# Normalize the view, otherwise the whole movies won't be fixed
# Output of the get_view command
cmd.set_view ("""\
     0.985269666,    0.099541061,    0.139044806,\
     0.137279570,    0.024403600,   -0.990231931,\
    -0.101961575,    0.994734287,    0.010378729,\
     0.000000000,    0.000000000, -258.863006592,\
    41.115470886,   41.043762207,    6.654131889,\
   204.089569092,  313.636444092,  -20.000000000 
""")


# Read trajectory in
trj = Trajectory("md.xtc")

# This requires having set up the GMX_DLL variable pointing to gromacs shared libraries

# Iterate over frames
for frame in trj:
    # It looks like we might need to write evrything to file, since the load_coords is not usable as is
    cmd.load(fn)

    cmd.ray()
    cmd.png('png/' + fn[:-4] + '.png')

    # Deleting the model should be quicker than relaunching pymol
    cmd.delete(fn[:-4])

# End Pymol Session
cmd.quit()


# Then make a movie with the pngs ...
# Source : https://trac.ffmpeg.org/wiki/Encode/H.264
# Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702



#DEBUG
#pymol.cmd.load_coords([ [0.8,0.8,0.8], [4.0, 0.0, 0.0] ], fn[:-4] , 1)
