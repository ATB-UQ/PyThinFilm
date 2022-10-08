import pmx
import math
import logging
from functools import partial
from PyThinFilm.helpers import atomic_density, insert_residues_faster, remove_residues_faster
import numpy as np

# Helper function for choosing which residues to insert
# Only accept residues with all atoms between min_z and max_z
def accept_residue(min_z, max_z, res):
    for a in res.atoms:
        if a.x[2] < min_z or a.x[2] > max_z:
            return False
    return True

# Handler class
class InsertionHandler(object):

    def __init__(self, model, split_min, split_max, density_thresh, bin_sz, split_ht=0, split_ht_tol=0):
        if isinstance(model, str):
            self.gro_file = model
            self.reload_model()
        else:
            self.gro_file = None
            self.model = model
        self.density_profile = None
        self.split_min = split_min
        self.split_max = split_max
        self._valid_residues_cache = None
        self.density_thresh = density_thresh
        self.bin_sz = bin_sz
        self.split_ht = split_ht
        self.split_ht_tol = split_ht_tol
        self.ranked_layers = None
        self._residues_translated = None

    def reload_model(self):
        logging.debug(f"Loading insertion handler model from file: {self.gro_file}")
        self.model = pmx.Model(self.gro_file)
        self.model.a2nm()

    def calc_density_profile(self, exclude_residues: list=[], set_profile=None):
        """
        Calculate and cache atomic density profile of solute.
        If `set_profile` is not `None`, density profile is set directly to that instead of calculating."
        """
        if set_profile is None:
            logging.debug("Calculating density profile")
            self.density_profile = atomic_density(self.model, self.bin_sz, exclude_residues)
        else:
            self.density_profile = set_profile

    def analyse_splits(self):
        """
        Search for consecutive bins between `self.split_min` and
        `self.split_max` with a solute density lower than the threshold and
        return a list of their z values (or an empty list if no valid splits
        are found).

        Requires that `calc_density_profile()` has been called.
        """
        if self.split_max <= self.split_min:
            logging.error("Invalid height range for splitting points")
            return []

        valid_splits = []
        assert(self.density_profile is not None)

        b = int(math.floor(self.split_min/self.bin_sz))
        b = 1 if b < 1 else b
        last_b = int(math.ceil(self.split_max/self.bin_sz))
        last_b = len(self.density_profile)-1 if last_b >= len(self.density_profile) else last_b
        logging.debug(f"Looking for splitting points between indices {b} and {last_b} of {len(self.density_profile)-1}")
        while b <= last_b:
            # Require two consecutive bins that contain only the specified
            # residue (or nothing at all).
            # Splitting point is between the two bins
            if not self.density_profile[b] < self.density_thresh:
                b += 2
            elif not self.density_profile[b-1] < self.density_thresh:
                b += 1
            else:
                valid_splits.append(b*self.bin_sz)
                b += 1
        return valid_splits

    def get_best_layers(self):
        """
        Find all combinations of splitting points that give the desired height within +/- tolerance,
        where the tolerance is `self.split_ht_tol * self.split_ht`.

        Output is a sorted list of tuples containing `(layer_bottom, layer_top, rating)`, where
        `rating` is the absolute difference between the layer height and `self.split_ht`
        """
        if self.ranked_layers is None:
            if self.split_max - self.split_min < self.split_ht * (1-self.split_ht_tol):
                logging.error(f"Not enough space between {self.split_min} and {self.split_max} to create a layer of height {self.split_ht}.")
                return []

            valid_splits = self.analyse_splits()
            if len(valid_splits) <= 1:
                logging.error(f"No possible layers found between {self.split_min} and {self.split_max}")
                return []
            if (max(valid_splits) - min(valid_splits)) < self.split_ht * (1-self.split_ht_tol):
                logging.error(f"No valid layers large enough. Maximum available is height is {max(valid_splits) - min(valid_splits)}")
                return []

            # Find combinations of values in valid_splits that give
            # a height within tolerance of split_ht
            ranked_layers = []
            last_min = 1
            split_zmax = max(valid_splits)
            for z_min in valid_splits:
                if (split_zmax - z_min) < self.split_ht * (1-self.split_ht_tol): break
                for i_zmax in range(last_min, len(valid_splits)):
                    z_max = valid_splits[i_zmax]
                    if (z_max - z_min) < self.split_ht * (1-self.split_ht_tol):
                        last_min = i_zmax
                        continue
                    if (z_max - z_min) > self.split_ht * (1+self.split_ht_tol):
                        break

                    ranked_layers.append( (z_min, z_max, abs((z_max-z_min) - self.split_ht)) )

            logging.debug(f"Found {len(ranked_layers)} possible layers for insertion")
            if len(ranked_layers) == 0:
                return []

            ranked_layers.sort(key=lambda x:x[2])
            self.ranked_layers = ranked_layers
        return self.ranked_layers

    def choose_random_layer(self, rng, strategy="weighted"):
        """Choose a random layer between two split points using the specified strategy."""
        if strategy == "best":
            return self.get_best_layers()[0]
        elif strategy == "random":
            layers = self.get_best_layers()
            return rng.sample(layers, 1)[0]
        else:
            if strategy != "weighted":
                logging.warning(f"Unknown layer insertion strategy \"{strategy}\". Defaulting to \"weighted\".")
            layers = self.get_best_layers()
            weights = np.array(layers)[:,2]
            max_wt = np.max(weights)
            if max_wt == 0:
                max_wt = 1
            weights = max_wt - weights
            weights = np.cumsum(weights)
            weights /= weights[-1]
            random_num = rng.random()
            layer_id = 0
            while weights[layer_id] < random_num:
                layer_id += 1
            return layers[layer_id]

    def get_valid_residues(self):
        """Return and cache residues between `self.split_min` and `self.split_max`"""
        if self._valid_residues_cache is None:
            # This loops over all residues so could be very slow.
            # Using map so that it can parallelise easily later if needed (e.g. with multiprocessing)
            # Deep copy residues in case using self
            self._valid_residues_cache = \
                [r for r, keep in zip(
                    self.model.residues,
                    map(partial(accept_residue,
                                self.split_min,
                                self.split_max),
                                self.model.residues)
                    ) if keep
                 ]
        return self._valid_residues_cache

    def insert(self, insert_z, aux_solution, layer_min, layer_max, substrate, extra_space = 0.1):
        """Insert all residues in `input_model` between `layer_min` and `layer_max` into a gap created above insert_z."""
        # Restore model in case this isn't the first insertion from it.
        aux_solution.restore_model()

        # Find residues in slab to be inserted. Do this first in case using self.
        # Loop over only cached residues to save time
        logging.debug("Finding valid residues for insertion")
        valid_residues = aux_solution.get_valid_residues()
        logging.debug("Filtering residues for insertion")
        res_to_add = \
            [r for r, keep in zip(
                valid_residues,
                map(partial(accept_residue,
                            layer_min,
                            layer_max),
                    valid_residues)
                ) if keep
            ]
        logging.debug(f"Found {len(res_to_add)} residues to add")

        # Make space in system
        space = layer_max - layer_min + 2*extra_space
        self.model.box[2][2] += space
        logging.debug(f"    Making {space} nm of space")
        res_to_remove = []
        for res in self.model.residues:
            if res.resname == substrate:
                # Move any substrate atoms clipping through z boundary to
                # account for change in box size
                for atom in res.atoms:
                    # This test could cause issues if insert_z is very close to
                    # the substrate, but that wouldn't be desired anyway since
                    # there are probably interesting things going on in that region
                    if atom.x[2] >= insert_z:
                        atom.translate([0, 0, space])
            else:
                # Move any residue that has an atom higher than the insertion
                # height, so extra_space can be small (~ 2xVDW radius)
                move_res = False
                for atom in res.atoms:
                    if atom.x[2] >= insert_z:
                        move_res = True
                        break
                if move_res:
                    # Remove any residues that cross the plane and translate those above it
                    remove = False
                    for atom in res.atoms:
                        if atom.x[2] < insert_z:
                            remove = True
                            res_to_remove.append(res)
                            break
                    if not remove:
                        res.translate([0, 0, space])

        delta_mixture = {}
        if len(res_to_remove) > 0:
            logging.debug(f"Removing residues that cross the plane z = {insert_z}")
            remove_residues_faster(self.model, res_to_remove)
            for res in res_to_remove:
                if res.resname not in delta_mixture:
                    delta_mixture[res.resname] = 0
                delta_mixture[res.resname] -= 1

        # Translate residues to the correct z
        delta_z =  insert_z - layer_min + extra_space
        for res in res_to_add:
            res.translate([0, 0, delta_z])
            if res.resname not in delta_mixture:
                delta_mixture[res.resname] = 0
            delta_mixture[res.resname] += 1

        # Register the change in the auxiliary system so we can undo it later.
        # Much faster to do this than to reload model from file.
        aux_solution.register_change(res_to_add, delta_z)

        # Insert the new residues
        insert_residues_faster(self.model, res_to_add)

        return delta_mixture

    def register_change(self, residues_translated, delta_z):
        """Register translation of molecules due to an insertion"""
        self._residues_translated = {
            "residues": residues_translated,
            "delta_z": delta_z
        }

    def restore_model(self):
        """Restore translations due to a previous insertion"""
        if self._residues_translated is None:
            return
        delta_z = -self._residues_translated["delta_z"]
        for res in self._residues_translated["residues"]:
            res.translate([0, 0, delta_z])
        self._residues_translated = None

