import pmx
import cPickle
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d

MODEL_PATH_TEMPLATE = "/mddata2/uqmstroe/OLED_data/{0}wpc_end.gro"
MODEL_CACHE_TEMPLATE = "cached_modles.pickle"
USE_CACHE = False
WPC = (2, )#6, 15, 30)

BOX_IMAGE_SIGNS = 
#print BOX_IMAGE_SIGNS

def get_box_image_vectors(box):
    return [
        vec([box_dim * sign for (sign, box_dim) in zip(signs, box)])
        for signs in BOX_IMAGE_SIGNS
    ]
def pbc_translations(x, y, z):
    images = list(product([1, 0, -1], repeat=3))
    return np.array([x,y,z]) * np.array(images)

def distance(a, b):
    return np.linalg.norm(np.array(a)-np.array(b))

def pbc_distances(points, box_x, box_y, box_z):
    distances = []
    for translation in pbc_translations(box_x, box_y, box_z):
        dists.append(np.sqrt((x+translation[0]-center[0])**2 + (y+translation[1]-center[1])**2))
    return any([dist <= RADIUS for dist in dists])

def load_model(model_path):
    model = pmx.Model(model_path)
    model.nm2a()
    return model

def load_models():
    if USE_CACHE and exists(MODEL_CACHE_TEMPLATE):
        print "Loading models from cache"
        with open(MODEL_CACHE_TEMPLATE) as fh:
            models = cPickle.load(fh)
        print "Loaded {0}-wpc modles".format(",".join(map(str, models.keys())))
        return models
    models = {}
    for wpc in WPC:
        print "Loading {0}wpc model".format(wpc)
        model_path = MODEL_PATH_TEMPLATE.format(wpc)
        models[wpc] = load_model(model_path)
    if USE_CACHE:
        with open(MODEL_CACHE_TEMPLATE, "w") as fh:
            cPickle.dump(models, fh, protocol=cPickle.HIGHEST_PROTOCOL)
    return models

def maxZHeight(model):
        return max([a.x[2] for a in model.atoms])

def molecule_number(model):
        return len(model.residues)

def filter_model(model, keep_molecules):
    keeping = []
    for molecule in model.residues:
        if molecule.resname in keep_molecules:
            keeping.append(molecule)
    model.residues = keeping

def get_coms(model):
    coms = np.empty((len(model.residues), 3), dtype=float)
    for i, molecule in enumerate(model.residues):
        coms[i] = molecule.com(vector_only=True)
    return coms

def scatter_3d(xs, ys, zs):

    #xmin = min(xs); xmax = max(xs)
    #ymin = min(ys); ymax = max(ys)

    fig = plt.figure()
    ax = mpl_toolkits.mplot3d.Axes3D(fig)
    ax.scatter(xs, ys, zs)
    plt.show()

def main():
    models = load_models()
    for model in models.values():
        print model
        filter_model(model, ["IPR", "IPS"])
        coms = get_coms(model)
        print coms[0]
        xs, ys, zs = zip(*coms)
        pbc_distance_matrix = pbc_distances(coms, model.box[0][0], model.box[1][1], model.box[2][2])
        scatter_3d(xs, ys, zs)
        for molecule in model.residues:
            for atom in molecule.atoms:
                print atom

if __name__=="__main__":
    main()