import logging
import pathlib
from pathlib import Path
import numpy as np
from math import erf

from PyThinFilm.common import RESOURCES_DIR


def foldnorm_mean(m, s):
    """https://en.wikipedia.org/wiki/Folded_normal_distribution"""
    return s*np.sqrt(2/np.pi)*np.exp((-m**2)/(2*s**2)) + m*erf(m/np.sqrt(2*s**2))


def res_highest_z(residue):
    return max(atom.x[2] for atom in residue.atoms)


def remove_residues_faster(model, residues):

    assert len(model.chains) == 1
    chain = model.chains[0]
    for residue in residues:
        logging.debug("Removing residue with id: {} ({})".format(residue.id, residue.resname))
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


def insert_residues_faster(model, residues):
    """Faster insertion of multiple residues"""
    assert len(model.chains) == 1
    chain = model.chains[0]
    for res in residues:
        chain.residues.append(res)
        chain.model.residues.append(res)
    model.atoms = []
    for r in model.residues:
        for a in r.atoms:
            model.atoms.append(a)
    chain.atoms = []
    for r in chain.residues:
        for a in r.atoms:
            chain.atoms.append(a)
    model.renumber_atoms()
    model.renumber_residues()
    chain.make_residue_tree()


def recursive_correct_paths(node):
    for key, value in node.items():
        if isinstance(value, dict):
            recursive_correct_paths(value)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    recursive_correct_paths(item)
        elif "_file" in key:
            orig_path = Path(value)
            if orig_path.exists():
                node[key] = orig_path.absolute()
            else:
                resource_path = RESOURCES_DIR/value
                if not resource_path.exists():
                    logging.warning(f"File not found: {key} = {value}")
                node[key] = resource_path
                logging.debug(f"Path updated: {key} = {resource_path}")


def recursive_convert_paths_to_str(node):
    for key, value in node.items():
        if isinstance(value, dict):
            recursive_convert_paths_to_str(value)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    recursive_convert_paths_to_str(item)
        elif isinstance(value, pathlib.Path):
            node[key] = value.as_posix()


def get_mass_dict(itp_string):
    itp_string = itp_string.split("[ atoms ]")[1].split("[ bonds ]")[0]
    mass_dict = {}
    for line in itp_string.splitlines():
        if not line or line.startswith(";"):
            continue
        mass_dict[line.split()[4]] = float(line.split()[7])
    return mass_dict


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


def mean_x(mol):
    """Calculate geometric center of a molecule."""
    return np.sum([a.x for a in mol.atoms], axis=0)/len(mol.atoms)


def mol_distance(mol1, mol2, box):
    """Calculate distance between geometric center of two molecules."""
    x1 = mean_x(mol1)
    x2 = mean_x(mol2)
    dx = [abs(p1 - p2) for p1, p2 in zip(x1, x2)]
    dx = [d if d < box[i][i]/2 else box[i][i] - d for i, d in enumerate(dx)]
    return np.sqrt(sum([d**2 for d in dx]))


def atomic_density(model, bin_sz, exclude_residues: list):
    """Calculate atomic density profile using bins of `bin_sz`, excluding atoms with resname in `exclude_residues`"""
    n_bins = int(np.ceil(model.box[2][2]/bin_sz))
    prof = np.zeros((n_bins+1,))
    for a in model.atoms:
        if a.resname not in exclude_residues:
            b = int(np.floor(a.x[2]/bin_sz))
            prof[b] += 1
    prof /= model.box[0][0]*model.box[1][1] * bin_sz
    return prof
