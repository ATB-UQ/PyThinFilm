import logging
from itertools import combinations
from pathlib import Path
import numpy as np
from math import erf

from PyThinFilm.common import RESOURCES_DIR


def foldnorm_mean(m, s):
    """https://en.wikipedia.org/wiki/Folded_normal_distribution"""
    return s*np.sqrt(2/np.pi)*np.exp((-m**2)/(2*s**2)) + m*erf(m/np.sqrt(2*s**2))


def res_highest_z(residue):
    return max(atom.x[2] for atom in residue.atoms)


def remove_residues(model,residues):
    logging.debug("num  residue: {0} num chains {1}".format(len(model.chains), len(model.residues)))
    assert len(model.chains) == 1
    chain = model.chains[0]
    for residue in residues:
        logging.debug("Removing residue: {0}".format(residue.resname))
        idx = chain.residues.index(residue)
        try:
            midx = chain.model.residues.index(residue)
        except:
            midx = -1
        del chain.residues[idx]
        del chain.model.residues[midx]
    model.atoms = []
    for r in model.residues:
        for atom in r.atoms:
            model.atoms.append(atom)
    chain.atoms = []
    for r in chain.residues:
        for atom in r.atoms:
            chain.atoms.append(atom)
    model.renumber_atoms()
    model.renumber_residues()
    chain.make_residue_tree()


def recursive_correct_paths(node):
    for key, value in node.items():
        if isinstance(value, dict):
            recursive_correct_paths(value)
        elif isinstance(value, list):
            for item in value:
                recursive_correct_paths(item)
        elif "_file" in key and not Path(value).exists():
            resource_path = RESOURCES_DIR/value
            if not resource_path.exists():
                logging.warning(f"File not found: {key} = {value}")
            node[key] = resource_path
            logging.debug(f"Path updated: {key} = {resource_path}")


def get_mass_dict(itpString):
    itpString = itpString.split("[ atoms ]")[1].split("[ bonds ]")[0]
    massDict = {}
    for line in itpString.splitlines():
        if not line or line.startswith(";"):
            continue
        massDict[line.split()[4]] = float(line.split()[7])
    return massDict


def group_residues(resnames):
    res_groups = []
    i = 0
    current_resname = resnames[0]
    if not resnames:
        return ""
    while current_resname != '':
        count = 1
        try:
            next_resname = resnames[i + 1]
        except IndexError:
            next_resname = ""
        while next_resname == current_resname :
            count +=1
            try:
                next_resname = resnames[i + count]
            except IndexError:
                next_resname = ""
        res_groups.append("{0} {1}".format(current_resname, count))
        i += count
        try:
            current_resname = resnames[i]
        except IndexError:
            current_resname = ""
    return res_groups


def basename_remove_suffix(p):
    return ".".join(Path(p).name.split('.')[0:2])


def calc_exclusions(model):
    Lx = model.box[0][0]
    Ly = model.box[1][1]
    exclusions = []
    atoms = [(a.id, np.array(a.x)) for residue in model.residues for a in residue.atoms]
    for ai, aj in combinations(atoms, 2):
        dx = ai[1][0] - aj[1][0]
        dy = ai[1][1] - aj[1][1]
        dx = dx - Lx if dx > 0.5 * Lx else dx
        dy = dy - Ly if dy > 0.5 * Ly else dy
        dx = dx + Lx if dx < -0.5 * Lx else dx
        dy = dy + Ly if dy < -0.5 * Ly else dy
        d = np.sqrt(dx * dx + dy * dy)
        if d < 0.3:
            # print(ai[0], aj[0], d)
            exclusions.append("{} {}".format(ai[0], aj[0]))
        # if len(exclusions) > 100:
        #     break

    with open("excl.txt", "w") as fh:
        fh.write("\n".join(exclusions))