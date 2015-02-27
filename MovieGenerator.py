import os
import shutil
import numpy
os.environ["GMX_DLL"]="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"  
import pmx
from pmx.xtc import Trajectory
from os.path import join, basename, exists
from subprocess import Popen, PIPE
import yaml
import argparse
import re
import glob
from copy import deepcopy
from jinja2 import Template
from operator import itemgetter

import logging
VERBOCITY = logging.DEBUG

from Deposition import GMX_PATH, TOPOLOGY_FILE

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

PYMOL_PNG_TEMPLATE="""
png {0}, width=1200, height=800, dpi=300, ray=1
"""

YAML_SCENES = glob.glob('scenes/*.yml')

FAST_PNG_LIMIT = 30

class MovieGenerator(object):

    def __init__(self, runConfig, runID):
        self.runConfig = runConfig
        self.dirname = runID
        self.sim_number = int(basename(runID))
        self.fn = self.absolute('init.gro')
        
        map(lambda x: os.makedirs(x) if not exists(x) else '', map(self.absolute, ['pml', 'pdb', 'png']))
        # Take the maximum of the frame_averaging over all the **active** scenes for a given sim_number
        self.average_n_frames = max( [ x["frame_averaging"] for x in map(lambda x: yaml.load(open(x)), YAML_SCENES) if x['first_sim_id'] <= self.sim_number <= x['last_sim_id'] ] )
        self.skip_n_frames = max( [ x["frame_skipping"] for x in map(lambda x: yaml.load(open(x)), YAML_SCENES) if x['first_sim_id'] <= self.sim_number <= x['last_sim_id'] ] )
        if self.skip_n_frames != 1 and self.average_n_frames != 1 :
            error_string = "Error: Can't have both frame averaging and frame skipping enabled at the same time. Set one of them to 1 and try again."
            logging.error(error_string)
            raise Exception(error_string)
        self.n_cores = self.runConfig['movies']['n_cores']
        
        self.fn_xtc_orig = self.absolute('md.xtc')
        self.fn_xtc = self.fn_xtc_orig if self.skip_n_frames == 1 else self.absolute('md_skipped_{0}.xtc'.format(self.skip_n_frames))  

    def absolute(self, path):
        return join(self.dirname, path)
    
    
    def fixPBCAndSkipFrames(self):
        args = "{GMX_PATH}trjconv_d -pbc mol -s md.tpr -f {xtc_in} -o temp_{xtc_out} -skip {skip} <<EOF\n0\nEOF".format(**{"GMX_PATH":GMX_PATH, 
                                                                                                                    "skip":self.skip_n_frames,
                                                                                                                    "xtc_out" : self.fn_xtc,
                                                                                                                    "xtc_in" : self.fn_xtc_orig,})
        logging.debug("running: {0}".format(args))
        Popen(args, shell=True, cwd=self.dirname).wait()
        
        shutil.move(join(self.dirname, "temp_{0}".format(self.fn_xtc)), join(self.dirname, self.fn_xtc))
        
    # pnx requires having set up the GMX_DLL variable pointing to gromacs shared libraries
    # GMX_DLL="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"
    
    # Iterate over frames

    def createPNG(self, m, tmp_base_fn, count):
        png_file = self.absolute("png/{0:0>4d}.png".format(count))
        tmp_fn = tmp_base_fn + str(count) + ".pdb"
        m.write(tmp_fn)
        pml_fn = self.absolute('pml/{0}.pml'.format(count))
        with open(pml_fn, 'w') as fp:
            # Scenes need to be read and written in numerial order
            fp.write( self.createPymolSceneString(tmp_fn) )
            # The PNG generation is the same for all of them
            fp.write(PYMOL_PNG_TEMPLATE.format(png_file))
        logging.debug("Processed png number {0} ({1})".format(count, png_file))
        Popen("pymol -qc -n {0}".format(pml_fn).split(), stdout=PIPE, stderr=PIPE).wait()
        #Popen("pymol -qc -n {0}".format(pml_fn).split()).wait()
        os.remove(tmp_fn)

    def createPymolSceneString(self, tmp_fn):
        strPML = ''
        # First, hard-code the number of threads pymol is allowed to use
        strPML += 'set max_threads, {0}\n'.format(self.n_cores)
        # Then, read and render YAML scene's pymol commands
        for sceneFile in YAML_SCENES :
            scene = yaml.load(open(sceneFile))
            # If the current scene is not in the scene range, continue to the next one
            if not (scene['first_sim_id'] <= self.sim_number <= scene['last_sim_id'] ) :
                continue
            # Reconstruct python string by joining the individual lines
            scenePML = '\n'.join(scene['pymol_commands'])
            # Interpret it with jinja2
            t = Template(scenePML)
            # Write it to file
            strPML += '\n' + t.render(tmp_fn=tmp_fn, last_resid=(self.sim_number+1) )
        return strPML

    def generatePNGs(self, fast_run=False):
        # First, flush all potential old png's
        self.flushPNGs()

        # Then, read trajectory in
        m = pmx.Model(self.fn)
        
        if self.average_n_frames > 1:
            mlist = [deepcopy(m) for _ in range(self.average_n_frames)]
        
        trj = Trajectory(self.fn_xtc)
        
        tmp_base_fn = self.absolute("pdb/" + basename(self.fn_xtc_orig).split(".")[0])
        
        count = 1
        if self.average_n_frames > 1:
            while self.shouldGenerateMorePNG(count, fast_run):
                # update each model with next from in trajectory    
                modelCounter = 0
                for _ in trj:
                    trj.update(mlist[modelCounter])
                    modelCounter += 1
                    if modelCounter == self.average_n_frames:
                        break
                # if there aren't self.average_n_frames number of models updated, exit while loop 
                if modelCounter != self.average_n_frames:
                    break
                
                if self.average_n_frames > 1:
                    for j, atom in enumerate(m.atoms):
                        atom.x[0] = numpy.mean([m_i.atoms[j].x[0] for m_i in mlist]) 
                        atom.x[1] = numpy.mean([m_i.atoms[j].x[1] for m_i in mlist])
                        atom.x[2] = numpy.mean([m_i.atoms[j].x[2] for m_i in mlist])
                self.createPNG(m, tmp_base_fn, count)
                count += 1
        else:
            for _ in trj:
                trj.update(m)
                self.createPNG(m, tmp_base_fn, count)
                count += 1
            
    def flushPNGs(self):
        # Remove all pngs (before regenerating them). This is useful when we modify the frame avaraging and end up overwriting less pngs than previously generated.
        pngDir = join(self.dirname,"png")
        if exists(pngDir) :
            logging.warning("Removing tree: {0}".format(pngDir))
            shutil.rmtree(pngDir)
            # Then, make the png directory back, just in case.
            os.makedirs(pngDir)
            
    def shouldGenerateMorePNG(self, count, fast_run):
        return (not fast_run) or (count <= FAST_PNG_LIMIT )

    # Then make a movie with the pngs ...
    # Source : http://robotics.usc.edu/~ampereir/wordpress/?p=702
    def generateMovie(self, fast_run=False, keep_png=False):
        # Log
        logging.info("Running movie generation for: {0}".format(self.dirname))

        abs_md_mp4 = self.absolute("md.mp4")
        # First, clean up potential old "md.mp4"
        # This is necessary since we will use the presence of this file as the proof ffmpeg succeeded
        if exists(abs_md_mp4) : os.remove(abs_md_mp4)
        # Then, generate the overlaid text
        overlaid_command = self.FFmpegOverlaidCommand()
        fast_run_command = "" #if not fast_run else "-r 1"
        args = "ffmpeg {3} -i {0} {2} {1}".format(self.absolute(r'png/%04d.png'), abs_md_mp4, overlaid_command, fast_run_command )
        logging.debug("running: {0}".format(args))
        p = Popen(args, shell=True, stdout=PIPE, stderr=PIPE)
        p.wait()
        error = p.stderr.read()
        # ffmpeg writes everything to stderr, so we have to deal with it manually by checking the file exists...
        if not exists( abs_md_mp4 ) : 
            logging.error("Generating movie from png with ffmpeg failed with error message: {0}. Check the log.".format(error))
        else :
            if not keep_png: self.flushPNGs()

    def FFmpegOverlaidCommand(self):
        counter_text = "[in]"
        BASE_HEIGHT = 25
        LINE_HEIGHT = 45
        height = BASE_HEIGHT - LINE_HEIGHT
        for compound, counter in self.getSortedMixtureFromTopologyFile():
            height += LINE_HEIGHT
            if counter_text != "[in]" : counter_text += ", "
            counter_text += "drawtext=fontfile=/home/uqbcaron/.fonts/OpenSans-Regular.ttf:text='{0} {1}':fontsize=40:x=25:y={2}".format(compound, counter, height)
        return '-vf "{0}"'.format(counter_text)

    def getSortedMixtureFromTopologyFile(self):
        mixture = {}
        with open( join(self.dirname, TOPOLOGY_FILE) ) as fh:
            for line in fh:
                m = re.search("^([a-zA-Z0-9]+) ([0-9]+)", line)
                if m :
                    mixture.setdefault(m.group(1), 0)
                    mixture[m.group(1)] += int(m.group(2))
        # Don't forget to remote substrate name
        del mixture[ self.runConfig["substrate"]["res_name"] ]
        # Return the list sorted by alphabetical order
        return sorted(mixture.items(), key=itemgetter(0))
    
def concatenateMovies(fileList, workdir):
    with open(join(workdir, "file.list"), 'w') as f:
        for filename in fileList:
            f.write("file {0}\n".format(filename))
            
    args = "yes | ffmpeg -f concat -i file.list md_tot.mp4"
    logging.debug("running: {0}".format(args))
    #Popen(args.split(), stdout=PIPE, stderr=PIPE).wait()
    Popen(args, shell=True, cwd=workdir).wait()

def isRunDir(f):
    return os.path.isdir(f) and re.search(r"^[0-9]+$", basename(f))

def mp4Exists(dirname):
    return "md.mp4" in os.listdir(dirname)

def getMovieFileList(workdir, batchStr):
    
    globPattern = "{0}/[0-9]*/md.mp4".format(workdir)
    logging.debug("globbing pattern to find movie files: " + globPattern)
    fileList = sorted(glob.glob(globPattern), key=lambda x:int(x.split("/")[-2]))
    
    if not batchStr: return fileList
    
    start, end = map(int, batchStr.split(":"))
    runRange = range(start, end + 1)
    
    filteredFileList = filter(lambda x:int(x.split("/")[-2]) in runRange, fileList)
    
    return filteredFileList

def getSortedRunDirList(dirname, batchStr):
    if not batchStr: return sorted(filter(isRunDir, map(lambda x:join(dirname, x), os.listdir(dirname))), key=lambda x:int(basename(x)))
    
    start, end = map(int, batchStr.split(":"))
    runRange = range(start, end+1)
    
    fullList = sorted(filter(isRunDir, map(lambda x:join(dirname, x), os.listdir(dirname))), key=lambda x:int(basename(x)))
    existsList = filter(lambda x: int(basename(x)) in runRange, fullList)
    
    if len(existsList) != len(runRange):
        logging.warning("Some runs specified in the range do not (yet) exist. Running those that do.")
                        
    return existsList

def generateAllIndividualMovies(runConfig, workdir, args):
    logging.basicConfig(level=VERBOCITY, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
   
    sortedRunDirList = getSortedRunDirList(workdir, args.batch)
    logging.info("Will run movie generation on the following sorted directory list :{0}".format(sortedRunDirList))
    for runID in sortedRunDirList :
        
        movie_generator = MovieGenerator(runConfig, runID)
        
        if mp4Exists(movie_generator.dirname) and not args.overwrite:
            continue
        
        if args.skip_png_generation:
            logging.warning("Skipping png generation in directory: {0}".format(runID))
        else:
            movie_generator.fixPBCAndSkipFrames()
            movie_generator.generatePNGs(fast_run=args.fast_run)
            
        movie_generator.generateMovie(fast_run=args.fast_run, keep_png=args.keep_png)

def parseCommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-d', '--overwrite', dest='overwrite', action='store_true')
    parser.add_argument('-b', '--batch', dest='batch')
    parser.add_argument('-c', '--concatenate', dest='concatenate', action='store_true')
    parser.add_argument('-sp', '--skip_png_generation', dest='skip_png_generation', action='store_true')
    parser.add_argument('-f', '--fast', dest='fast_run', action='store_true')
    parser.add_argument('-kp', '--keep_png', dest='keep_png', action='store_true')
    
    
    args = parser.parse_args()
    
    runConfig = yaml.load(open(args.input))
    workdir = join(PROJECT_DIR, runConfig["work_directory"])
    
    generateAllIndividualMovies(runConfig, workdir, args)
    
    if args.concatenate:
        logging.info("Launching sub-movie concatenation in: {0}".format(workdir))
        concatenateMovies(getMovieFileList(workdir,args.batch), workdir)

if __name__=="__main__":
    parseCommandline()
    
