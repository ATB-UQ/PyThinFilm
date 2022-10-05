import pmx
import yaml
import argparse
import math
import random
import logging
from functools import partial
from PyThinFilm.helpers import insert_residues_faster, remove_residues_faster

# Definition here required for parallel residue checking
def accept_residue(min_h, max_h, res):
    for a in res.atoms:
        if a.x[2] < min_h or a.x[2] > max_h: return False
    return True

# Handler class
class InsertionHandler(object):

    def __init__(self, model, layer_height = 0.5, from_model = False):
        if from_model:
            self.model = model
        else:
            self.model = pmx.Model(model)
            self.model.nm2a()
        self.counts = []
        self.layer_height = layer_height
        self.get_counts()

    # Find number of atoms from each residue type in each layer
    def get_counts(self):
        numBins = int(math.ceil(self.model.box[2][2]*10/self.layer_height))
        self.counts = list({} for _ in range(numBins+1))
        for a in self.model.atoms:
            b = int(math.floor(a.x[2]/self.layer_height))
            if a.resname not in self.counts[b]:
                self.counts[b][a.resname] = 0
            self.counts[b][a.resname] += 1

    # Helper function to check for layers that contain
    # < thresh of undesired residue atoms
    def is_valid_split(self, b, resname, thresh):
        if len(self.counts[b]) == 0:
            return True
        elif resname not in self.counts[b]:
            return False
        c = 0
        for key in self.counts[b]:
            if key != resname:
                c += self.counts[b][key]
        if float(c)/(c+self.counts[b][resname]) > thresh:
            #logging.debug(f"{100*float(c)/(c+self.counts[b][resName])}% non solvent atoms > {thresh}")
            return False
        else:
            return True

    # Search for consecutive layers that only include resname
    # and return a list of z values at the midpoint between 
    # those layers
    def analyse_splits(self, resname, min_h, max_h, thresh):
        if max_h <= min_h:
            msg = "Invalid height range"
            logging.error(msg)
            raise AssertionError(msg)
        validSplits = []

        #validLayer = lambda b: len(self.counts[b]) == 0 or (len(self.counts[b])==1 and self.counts[b].has_key(resName))

        b = int(math.floor(min_h/self.layer_height))
        b = 1 if b < 1 else b
        last_b = int(math.ceil(max_h/self.layer_height))
        last_b = len(self.counts)-1 if last_b >= len(self.counts) else last_b
        logging.debug(f"Searching from indices {b} to {last_b} of {len(self.counts)-1}")
        while b <= last_b:
            # Require two consecutive layers that contain only the specified residue (or nothing at all)
            if not self.is_valid_split(b, resname, thresh):
                b += 2
            elif not self.is_valid_split(b-1, resname, thresh):
                b += 1
            else:
                validSplits.append(b*self.layer_height)
                b += 1

        return validSplits


    # Find all splits between minH and maxH that gives the desired height within +/- tolerance
    def get_best_splits(self, resname, min_h, max_h, split_h, tol, thresh):
        assert (max_h - min_h) >= split_h * (1-tol), f"Not enough space between {min_h} and {max_h} for a layer of height {split_h}"
        valid_splits = self.analyse_splits(resname, min_h, max_h, thresh)

        if len(valid_splits) <= 1:
            msg = f"No possible splits found between {min_h} and {max_h}"
            logging.error(msg)
            raise AssertionError(msg)
        if (max(valid_splits) - min(valid_splits)) < split_h * (1-tol):
            msg = f"No valid splits large enough! Maximum available is height is {max(valid_splits) - min(valid_splits)}"
            logging.error(msg)
            raise AssertionError(msg)

        ranked_splits = []
        last_min = 1
        for z_min in valid_splits:
            if (max(valid_splits) - z_min) < split_h * (1-tol): break
            for i_zMax in range(last_min, len(valid_splits)):
                z_max = valid_splits[i_zMax]
                if (z_max - z_min) < split_h * (1-tol):
                    last_min = i_zMax
                    continue
                if (z_max - z_min) > split_h * (1+tol): break

                ranked_splits.append([z_min, z_max, abs((z_max-z_min) - split_h)])

        assert len(ranked_splits) > 0, "No valid splits found"
        logging.info(f"Found {len(ranked_splits)} possible splits")

        # TODO: rank by closest to desired solute ratio?
        ranked_splits.sort(key=lambda x:x[2])

        return ranked_splits


    # Choose either a random split or the best split
    def choose_split(self, ranked_splits, best = True):
        if best:
            return ranked_splits[0]
        else:
            return random.sample(ranked_splits, k=1)[0]

    # Insert all residues in model between minH and maxH 
    # into a gap created at insertH
    def insert(self, insert_h, input_model, min_h, max_h, substrate, extra_space = 0.1): #, topSpace = 1):
        # Make space in system
        self.model.box[2][2] += (max_h - min_h + 2*extra_space) * 0.1
        logging.debug(f"    Making {max_h - min_h + 2*extra_space} A of space")
        res_to_remove = []
        for res in self.model.residues:
            # Move any substrate atoms clipping through z boundary to account for change in box size
            if res.resname == substrate:
                for atom in res.atoms:
                    if atom.x[2] >= insert_h:
                        atom.translate([0, 0, max_h - min_h + 2*extra_space])
            else:
                # Move any residue that has an atom higher than the insertion height, so extraSpace can be small (~ 2xVDW radius)
                moveRes = False
                for atom in res.atoms:
                    if atom.x[2] >= insert_h:
                        moveRes = True
                        break
                if moveRes:
                    res.translate([0, 0, max_h - min_h + 2*extra_space])
                    remove = False
                    for atom in res.atoms:
                        if atom.x[2] < insert_h + max_h - min_h + 2*extra_space:
                            remove = True
                            break
                    if remove:
                        res_to_remove.append(res)

        delta_mixture = {}
        # Remove residues that collide with inserted layer
        if len(res_to_remove) > 0:
            remove_residues_faster(self.model, res_to_remove)
            for res in res_to_remove:
                if res.resname not in delta_mixture:
                    delta_mixture[res.resname] = 0
                delta_mixture[res.resname] -= 1


        # Find residues in range
        logging.debug("    Finding residues and inserting")
        total = len(input_model.residues)
        logging.debug(f"        {total} residues to check")
        res_to_add = [r for r, keep in zip(input_model.residues, map(partial(accept_residue, min_h, max_h), input_model.residues)) if keep]
        logging.debug(f"        {len(res_to_add)} residues to add")


        # Add residues to model
        for res in res_to_add:
            res.translate([0, 0, insert_h - min_h + extra_space])
            if res.resname not in delta_mixture:
                delta_mixture[res.resname] = 0
            delta_mixture[res.resname] += 1
        insert_residues_faster(self.model, res_to_add)
        return delta_mixture

def run_soln_insertion(config, model=None, write_output=True):
    logging.debug(f"Sourcing from {config['input_min']} to {config['input_max']} A")
    logging.debug(f"Inserting between {config['insert_min']} to {config['insert_max']} A")

    # Find slab to insert
    source_GRO = InsertionHandler(config["input_gro"], layer_height=config["search_bin_size"])
    chosen_split = source_GRO.choose_split(
        source_GRO.get_best_splits(config["solvent"],
                               config["input_min"],
                               config["input_max"],
                               config["insert_thickness"],
                               config["thickness_tol"],
                               config["max_non_solvent_atoms"]
                               ),
        best=True if not "randomise" in config else not config["randomise"])
    logging.info(f"Taking layer from {chosen_split[0]} to {chosen_split[1]} A of source")

    # Find insertion height
    main_GRO = InsertionHandler(model if model is not None else config["system_gro"], layer_height=config["search_bin_size"], from_model=(model != None))
    valid_splits = main_GRO.analyse_splits(config["solvent"], config["insert_min"], config["insert_max"], config["max_non_solvent_atoms"])
    if len(valid_splits) == 0:
        msg = "No valid insertion points found in system! Try using a smaller bin size"
        logging.error(msg)
        raise AssertionError(msg)
    insert_h = random.sample(valid_splits, 1)[0]
    logging.info(f"Inserting at {insert_h} A with extra {config['extra_space']} A above and below")

    # Insert
    if not config.has_key("extra_space"): config["extra_space"] = 1.5
    delta_mixture = main_GRO.insert(insert_h, source_GRO.model, chosen_split[0], chosen_split[1], config["substrate"], config["extra_space"])

    # Finalise
    if write_output:
        logging.debug(f"Writing output: {config['output_gro']}")
        main_GRO.model.write(config["output_gro"], title=main_GRO.model.title)
        logging.debug("Done")

    return delta_mixture


def parse_cmd():
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

    run_soln_insertion(config)
if __name__ == "__main__":
    parse_cmd()
