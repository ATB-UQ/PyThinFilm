import pymol
import pymol.experimenting
from pmx.xtc import *
from pmx import *
from os.path import join, basename, dirname, abspath
import os
import time
import subprocess
import yaml

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

PYMOL_TEMPLATE="""\
set_view (\
     0.449250728,    0.289243400,   -0.845287681,\
    -0.884793460,    0.275103867,   -0.376111388,\
     0.123753794,    0.916872859,    0.379511178,\
    -0.000094324,    0.000004143, -346.935516357,\
    49.073009491,   49.164260864,   27.948474884,\
   188.773605347,  505.097320557,  -20.000000000 )
color black, resi 1
bg_color white
load {0}
png {1}, width=1200, height=800, dpi=300, ray=1
"""

class MovieGenerator(object):

  def __init__(self, runConfigFile, sim_number):
    self.runConfig = yaml.load(open(runConfigFile))
    self.dirname = join(PROJECT_DIR, self.run_config, sim_number)
    self.fn = absolute('init.pdb')
    self.fn_xtc = absolute('md.xtc')
    map(lambda x: os.makedirs(x) if not os.path.exists(x) else '', map(absolute, ['pml', 'pdb', 'png']))

  def absolute(path):
    return join(self.dirname, path)

  # pnx requires having set up the GMX_DLL variable pointing to gromacs shared libraries
  # GMX_DLL="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"

  # Iterate over frames
  def generatePNGs():
    # Read trajectory in
    m = Model(self.fn)
    trj = Trajectory(self.fn_xtc)

    tmp_base_fn = absolute("pdb/" + self.fn_xtc[:-4])
    for i, frame in enumerate(trj):
        trj.update(m)
        tmp_fn = absolute(tmp_base_fn + str(i) + ".pdb")
        m.write(tmp_fn)
        with open(absolute('pml/{0}.pml'.format(i)), 'w') as fp:
          fp.write(PYMOL_TEMPLATE.format(tmp_fn, absolute("png/{0:0>4d}.png".format(i))))
        subprocess.Popen("pymol -qc -n pml/{0}.pml".format(i).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        print "Processed png number {}".format(i)

  # Then make a movie with the pngs ...
  # Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702
  def generateMovie():
    subprocess.Popen("ffmpeg -i {}/%04d.png md.mp4".format(absolute('png')).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

def runMovieGenerator(runConfigFile):

  for i in map(lambda x:x+1,range(0,1)):
    movie_generator = MovieGenerator(runConfigFile, i)
    movie_generator.generatePNGs()
    movie_generator.generateMovie()

if __name__=="__main__":
    runMovieGenerator(sys.argv[1])
