import pmx
import yaml
import argparse
import math
import random
import logging
from functools import partial
from itertools import imap
from PyThinFilm.helpers import insert_residues_faster, remove_residues_faster

# Definition here required for parallel residue checking
def acceptResidue(minH, maxH, res):
    for a in res.atoms:
        if a.x[2] < minH or a.x[2] > maxH: return False
    return True

# Handler class
class groHandler(object):

    def __init__(self, model, layerHeight = 0.5, fromModel = False):
        if fromModel:
            self.model = model
        else:
            self.model = pmx.Model(model)
            self.model.nm2a()
        self.counts = []
        self.layerHeight = layerHeight
        self.getCounts()

    # Find number of atoms from each residue type in each layer
    def getCounts(self):
        numBins = int(math.ceil(self.model.box[2][2]*10/self.layerHeight))
        self.counts = list({} for b in range(numBins+1))
        for a in self.model.atoms:
            b = int(math.floor(a.x[2]/self.layerHeight))
            if not self.counts[b].has_key(a.resname):
                self.counts[b][a.resname] = 0
            self.counts[b][a.resname] += 1

    # Helper function to check for layers that contain
    # < thresh of undesired residue atoms
    def isValidSplit(self, b, resName, thresh):
        if len(self.counts[b]) == 0:
            return True
        elif not self.counts[b].has_key(resName):
            return False
        c = 0
        for key in self.counts[b]:
            if key != resName:
                c += self.counts[b][key]
        if float(c)/(c+self.counts[b][resName]) > thresh:
            #logging.debug("{0}% non solvent atoms > {1}".format(100*float(c)/(c+self.counts[b][resName]), thresh))
            return False
        else:
            return True

    # Search for consecutive layers that only include resName
    # and return a list of z values at the midpoint between 
    # those layers
    def analyseSplits(self, resName, minH, maxH, thresh):
        if maxH <= minH:
            msg = "Invalid height range"
            logging.error(msg)
            raise AssertionError(msg)
        validSplits = []

        #validLayer = lambda b: len(self.counts[b]) == 0 or (len(self.counts[b])==1 and self.counts[b].has_key(resName))

        b = int(math.floor(minH/self.layerHeight))
        b = 1 if b < 1 else b
        lastB = int(math.ceil(maxH/self.layerHeight))
        lastB = len(self.counts)-1 if lastB >= len(self.counts) else lastB
        logging.debug("Searching from indices {0} to {1} of {2}".format(b, lastB, len(self.counts)-1))
        while b <= lastB:
            # Require two consecutive layers that contain only the specified residue (or nothing at all)
            if not self.isValidSplit(b, resName, thresh): 
                b += 2
            elif not self.isValidSplit(b-1, resName, thresh):
                b += 1
            else:
                validSplits.append(b*self.layerHeight)
                b += 1

        return validSplits


    # Find all splits between minH and maxH that gives the desired height within +/- tolerance
    def getBestSplit(self, resName, minH, maxH, splitH, tol, thresh):
        assert (maxH - minH) >= splitH * (1-tol), "Not enough space between {0} and {1} for a layer of height {2}".format(minH, maxH, splitH)
        validSplits = self.analyseSplits(resName, minH, maxH, thresh)

        if len(validSplits) <= 1:
            msg = "No possible splits found between {0} and {1}".format(minH, maxH)
            logging.error(msg)
            raise AssertionError(msg)
        if (max(validSplits) - min(validSplits)) < splitH * (1-tol):
            msg = "No valid splits large enough! Maximum available is height is {0}".format(max(validSplits) - min(validSplits))
            logging.error(msg)
            raise AssertionError(msg)

        rankedSplits = []
        lastMin = 1
        for zMin in validSplits:
            if (max(validSplits) - zMin) < splitH * (1-tol): break
            for i_zMax in range(lastMin, len(validSplits)):
                zMax = validSplits[i_zMax]
                if (zMax - zMin) < splitH * (1-tol):
                    lastMin = i_zMax
                    continue
                if (zMax - zMin) > splitH * (1+tol): break

                rankedSplits.append([zMin, zMax, abs((zMax-zMin) - splitH)])

        assert len(rankedSplits) > 0, "No valid splits found"
        logging.info("Found {0} possible splits".format(len(rankedSplits)))

        # TODO: rank by closest to desired solute ratio?
        rankedSplits.sort(key=lambda x:x[2])

        return rankedSplits


    # Choose either a random split or the best split
    def chooseSplit(self, rankedSplits, best = True):
        if best:
            return rankedSplits[0]
        else:
            return random.sample(rankedSplits, k=1)[0]

    # Insert all residues in model between minH and maxH 
    # into a gap created at insertH
    def insert(self, insertH, inputModel, minH, maxH, substrate, extraSpace = 0.1): #, topSpace = 1):
        # Make space in system
        self.model.box[2][2] += (maxH - minH + 2*extraSpace) * 0.1
        logging.debug("    Making {0} A of space".format(maxH - minH + 2*extraSpace))
        resToRemove = []
        for res in self.model.residues:
            # Move any substrate atoms clipping through z boundary to account for change in box size
            if res.resname == substrate:
                for atom in res.atoms:
                    if atom.x[2] >= insertH:
                        atom.translate([0, 0, maxH - minH + 2*extraSpace])
            else:
                # Move any residue that has an atom higher than the insertion height, so extraSpace can be small (~ 2xVDW radius)
                moveRes = False
                for atom in res.atoms:
                    if atom.x[2] >= insertH:
                        moveRes = True
                        break
                if moveRes:
                    res.translate([0, 0, maxH - minH + 2*extraSpace])
                    remove = False
                    for atom in res.atoms:
                        if atom.x[2] < insertH + maxH - minH + 2*extraSpace:
                            remove = True
                            break
                    if remove:
                        resToRemove.append(res)

        deltaMixture = {}
        # Remove residues that collide with inserted layer
        if len(resToRemove) > 0:
            remove_residues_faster(self.model, resToRemove)
            for res in resToRemove:
                if not deltaMixture.has_key(res.resname):
                    deltaMixture[res.resname] = 0
                deltaMixture[res.resname] -= 1


        # Find residues in range
        logging.debug("    Finding residues and inserting")
        total = len(inputModel.residues)
        logging.debug("        {0} residues to check".format(total))
        resToAdd = [r for r, keep in zip(inputModel.residues, imap(partial(acceptResidue, minH, maxH), inputModel.residues)) if keep]
        logging.debug("        {0} residues to add".format(len(resToAdd)))


        # Add residues to model
        for res in resToAdd:
            res.translate([0, 0, insertH - minH + extraSpace])
            if not deltaMixture.has_key(res.resname):
                deltaMixture[res.resname] = 0
            deltaMixture[res.resname] += 1
        insert_residues_faster(self.model, resToAdd)
        return deltaMixture

def runSolnInsertion(config, model=None, writeOutput=True):
    logging.debug("Sourcing from {0} to {1} A".format(config["input_min"], config["input_max"]))
    logging.debug("Inserting between {0} to {1} A".format(config["insert_min"], config["insert_max"]))

    # Find slab to insert
    sourceGRO = groHandler(config["input_gro"], layerHeight=config["search_bin_size"])
    chosenSplit = sourceGRO.chooseSplit(
        sourceGRO.getBestSplit(config["solvent"],
                               config["input_min"],
                               config["input_max"],
                               config["insert_thickness"],
                               config["thickness_tol"],
                               config["max_non_solvent_atoms"]
                               ),
        best=True if not "randomise" in config else not config["randomise"])
    logging.info("Taking layer from {0} to {1} A of source".format(chosenSplit[0], chosenSplit[1]))

    # Find insertion height
    mainGRO = groHandler(model if model != None else config["system_gro"], layerHeight=config["search_bin_size"], fromModel=(model != None))
    validSplits = mainGRO.analyseSplits(config["solvent"], config["insert_min"], config["insert_max"], config["max_non_solvent_atoms"])
    if len(validSplits) == 0:
        msg = "No valid insertion points found in system! Try using a smaller bin size"
        logging.error(msg)
        raise AssertionError(msg)
    insertH = random.sample(validSplits, 1)[0]
    logging.info("Inserting at {0} A with extra {1} A above and below".format(insertH, config["extra_space"]))

    # Insert
    if not config.has_key("extra_space"): config["extra_space"] = 1.5
    deltaMixture = mainGRO.insert(insertH, sourceGRO.model, chosenSplit[0], chosenSplit[1], config["substrate"], config["extra_space"])

    # Finalise
    if writeOutput:
        logging.debug("Writing output: {0}".format(config["output_gro"]))
        mainGRO.model.write(config["output_gro"], title=mainGRO.model.title)
        logging.debug("Done")

    return deltaMixture


def parseCmd():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='YAML input configuration file')
    parser.add_argument('-s', '--source', help='Source .gro file to insert slab into')
    parser.add_argument('-o', '--output', help='.gro file to save output to.')
    parser.add_argument('-l', '--log', help='Log file')
    args = parser.parse_args()

    logging.basicConfig(filename=args.log, level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%d-%m-%Y %H:%M:%S')

    ymlInput = yaml.safe_load(open(args.input))
    config = ymlInput["solution_layer_insert"]  # Containing in one node to use as input from deposition simulation in future

    config['system_gro'] = args.source
    config['output_gro'] = args.output
    if config.has_key('use_self') and config['use_self'] == True:
        config['input_gro'] = args.source

    runSolnInsertion(config)
if __name__ == "__main__":
    parseCmd()
