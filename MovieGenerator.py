import os
import sys
os.environ["GMX_DLL"]="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"  
import pmx
from pmx.xtc import Trajectory
from os.path import join, basename
import subprocess
import yaml
import argparse
import re
import glob

import logging
VERBOCITY = logging.DEBUG

from Deposition import GMX_PATH

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

PYMOL_TEMPLATE="""
load {0}
set_view (\
     0.449250728,    0.289243400,   -0.845287681,\
    -0.884793460,    0.275103867,   -0.376111388,\
     0.123753794,    0.916872859,    0.379511178,\
    -0.000094324,    0.000004143, -346.935516357,\
    49.073009491,   49.164260864,   27.948474884,\
   188.773605347,  505.097320557,  -20.000000000 )
color black, resi 1
bg_color white
png {1}, width=1200, height=800, dpi=300, ray=1
"""

class MovieGenerator(object):

    def __init__(self, runConfig, sim_number):
        logging.info("Running movie generation for: {0}".format(sim_number))
        self.runConfig = runConfig
        self.dirname = join(PROJECT_DIR, self.runConfig["work_directory"], str(sim_number))
        self.fn = self.absolute('init.gro')
        self.fn_xtc = self.absolute('md_fixedPBC.xtc')
        map(lambda x: os.makedirs(x) if not os.path.exists(x) else '', map(self.absolute, ['pml', 'pdb', 'png']))

    def absolute(self, path):
        return join(self.dirname, path)
    
    
    def fixPBC(self):
        args = "{GMX_PATH}trjconv_d -pbc mol -s md.tpr -f md.xtc -o md_fixedPBC.xtc <<EOF\n0\nEOF".format(**{"GMX_PATH":GMX_PATH})
        logging.debug("running: {0}".format(args))
        subprocess.Popen(args, shell=True, cwd=self.dirname).wait()
        
    # pnx requires having set up the GMX_DLL variable pointing to gromacs shared libraries
    # GMX_DLL="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"
    
    # Iterate over frames
    def generatePNGs(self):
        # Read trajectory in
        m = pmx.Model(self.fn)
        trj = Trajectory(self.fn_xtc)
        
        tmp_base_fn = self.absolute("pdb/" + basename(self.fn_xtc)[:-4])
        for i, _ in enumerate(trj):
            png_file = self.absolute("png/{0:0>4d}.png".format(i))
            trj.update(m)
            tmp_fn = tmp_base_fn + str(i) + ".pdb"
            m.write(tmp_fn)
            pml_fn = self.absolute('pml/{0}.pml'.format(i))
            with open(pml_fn, 'w') as fp:
                fp.write(PYMOL_TEMPLATE.format(tmp_fn, png_file))
                
            logging.debug("Processed png number {0} ({1})".format(i, png_file))
            subprocess.Popen("pymol -qc -n {0}".format(pml_fn).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
            #subprocess.Popen("pymol -qc -n {0}".format(pml_fn).split()).wait()
            
    
    
    # Then make a movie with the pngs ...
    # Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702
    def generateMovie(self):
        args = "yes | ffmpeg -i {0} {1}".format(*[self.absolute(x) for x in (r'png/%04d.png','md.mp4')])
        logging.debug("running: {0}".format(args))
        #subprocess.Popen(args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        subprocess.Popen(args, shell=True).wait()
        
def concatenateMovies(fileList):
    with open("file.list", 'w') as f:
        for filename in fileList:
            f.write("file {0}\n".format(filename))
            
    args = "yes | ffmpeg -f concat -i file.list md_tot.mp4"
    logging.debug("running: {0}".format(args))
    #subprocess.Popen(args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
    subprocess.Popen(args, shell=True).wait()

def isRunDir(f):
    return os.path.isdir(f) and re.search(r"^[0-9]+$", basename(f))

def mp4Exists(dirname):
    return "md.mp4" in os.listdir(dirname)

def getMovieFileList(workdir):
    globPattern = "{0}/[0-9]*/md.mp4".format(workdir)
    logging.debug("globbing pattern to find movie files: " + globPattern)
    fileList = glob.glob(globPattern)
    return fileList

def getSortedRunDirList(dirname):
    return sorted(filter(isRunDir, map(lambda x:join(dirname, x), os.listdir(dirname))), key=lambda x:int(basename(x)))

def runMovieGenerator(args):
    logging.basicConfig(level=VERBOCITY, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
    runConfig = yaml.load(open(args.input))
    workdir = join(PROJECT_DIR, runConfig["work_directory"])
    
    for runID in getSortedRunDirList(workdir):
        movie_generator = MovieGenerator(runConfig, runID)
        
        if mp4Exists(movie_generator.dirname) and not args.overwrite:
            continue
        
        movie_generator.fixPBC()
        movie_generator.generatePNGs()
        movie_generator.generateMovie()
    
    concatenateMovies(getMovieFileList(workdir))

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-D', '--overwrite', dest='overwrite', action='store_true')
    args = parser.parse_args()
    runMovieGenerator(args)
