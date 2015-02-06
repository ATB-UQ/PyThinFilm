import os
os.environ["GMX_DLL"]="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"  
import pmx
from pmx.xtc import Trajectory
from os.path import join, basename
import subprocess
import yaml
import sys

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

    def __init__(self, runConfig, sim_number):
        self.runConfig = runConfig
        self.dirname = join(PROJECT_DIR, self.runConfig["work_directory"], str(sim_number))
        self.fn = self.absolute('init.gro')
        self.fn_xtc = self.absolute('md.xtc')
        map(lambda x: os.makedirs(x) if not os.path.exists(x) else '', map(self.absolute, ['pml', 'pdb', 'png']))

    def absolute(self, path):
        return join(self.dirname, path)

    # pnx requires having set up the GMX_DLL variable pointing to gromacs shared libraries
    # GMX_DLL="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"
    
    # Iterate over frames
    def generatePNGs(self):
        # Read trajectory in
        m = pmx.Model(self.fn)
        trj = Trajectory(self.fn_xtc)
        
        tmp_base_fn = self.absolute("pdb/" + basename(self.fn_xtc)[:-4])
        for i, _ in enumerate(trj):
            trj.update(m)
            tmp_fn = tmp_base_fn + str(i) + ".pdb"
            m.write(tmp_fn)
            pml_fn = self.absolute('pml/{0}.pml'.format(i))
            with open(pml_fn, 'w') as fp:
                fp.write(PYMOL_TEMPLATE.format(tmp_fn, self.absolute("png/{0:0>4d}.png".format(i))))
            
            subprocess.Popen("pymol -qc -n {0}".format(pml_fn).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
            print "Processed png number {0}".format(i)
    
    
    # Then make a movie with the pngs ...
    # Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702
    def generateMovie(self):
        args = "ffmpeg -i {0} {1}".format(*[self.absolute(x) for x in (r'png/%04d.png','md.mp4')])
        print args
        subprocess.Popen(args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
def concatenateMovies(fileList):
    with open("file.list", 'w') as f:
        for filename in fileList:
            f.write("file {0}\n".format(filename))
            
    args = "ffmpeg -f concat -i file.list md_tot.mp4"
    print args
    subprocess.Popen(args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

def runMovieGenerator(runConfigFile):
    runConfig = yaml.load(open(runConfigFile))
    fileList = []
    for runID in (1,2):
        movie_generator = MovieGenerator(runConfig, runID)
        fileList.append(movie_generator.absolute('md.mp4'))
        movie_generator.generatePNGs()
        movie_generator.generateMovie()
    
    concatenateMovies(fileList)

if __name__=="__main__":
    runMovieGenerator(sys.argv[1])
