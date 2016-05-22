import cPickle
from os.path import exists
import numpy as np
from numpy import histogram
import matplotlib.pyplot as plt
from itertools import product, combinations
import mpl_toolkits.mplot3d
from scipy.spatial.distance import pdist, squareform
from graph_tool.all import Graph, graph_draw
from graph_tool.search import bfs_iterator
from graph_tool.stats import vertex_average, vertex_hist
import pmx

MODEL_PATH_TEMPLATE = "{0}wpc_end.gro"
MODEL_CACHE_TEMPLATE = "cached_modles.pickle"
USE_CACHE = False
WPC = (2, 6, 15, 30)

DIMENSIONS = 3
N_IMAGES = DIMENSIONS**2

def pbc_translations(x, y, z):
    # generate all unity periodic image translations, sorted such that the identity transformation [0,0,0] is first
    images = sorted(product([1, 0, -1], repeat=DIMENSIONS), key=lambda x:np.sum(np.abs(x)))
    return np.array([x,y,z]) * np.array(images)

def xy_pbc_translations(x, y, z):
    images = sorted([t for t in product([1, 0, -1], repeat=DIMENSIONS) if t[2] == 0], key=lambda x:np.sum(np.abs(x)))
    return np.array([x,y,z]) * np.array(images)

def distance(a, b):
    return np.linalg.norm(np.array(a)-np.array(b))

def get_pbc_points(points, box_x, box_y, box_z):
    pbc_points = np.empty((N_IMAGES, points.shape[0], DIMENSIONS))
    for i, translation in enumerate(xy_pbc_translations(box_x, box_y, box_z)):
        pbc_points[i] = points + translation
    return np.reshape(pbc_points, (points.shape[0]*N_IMAGES, DIMENSIONS)) # flatten to 1D array of (x,y,z) elements

def get_explicit_pbc_distances(points, box):
    box_x, box_y, box_z = box[0][0], box[1][1], box[2][2]
    pbc_points = get_pbc_points(points, box_x, box_y, box_z)
    distances = pdist(pbc_points)
    distances = squareform(distances)
    return pbc_points, distances

def get_implicit_pbc_distances(points, box):
    distances = np.empty((points.shape[0], points.shape[0]))
    for i,j in product(range(points.shape[0]), repeat=2):
        distances[i][j] = min_pbc_distance(points[i], points[j], box)
    return distances

def min_pbc_distance(point1, point2, box):
    box_x, box_y, box_z = box[0][0], box[1][1], box[2][2]
    pbc_points = get_pbc_points(np.array([point1, point2]), box_x, box_y, box_z)
    return np.min(pdist(pbc_points))

def load_model(model_path):
    model = pmx.Model(model_path)
    #model.nm2a()
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

def scatter_3d(xs, ys, zs, box=None, distance_cutoff_graph=None, cutoff_distance=None):

    #xmin = min(xs); xmax = max(xs)
    #ymin = min(ys); ymax = max(ys)

    fig = plt.figure()
    ax = mpl_toolkits.mplot3d.Axes3D(fig)
    ax.scatter(xs, ys, zs, alpha=0.3)
    if box:
        box = [box[0][0], box[1][1], box[2][2]]
        unit_box_corner_vectors = list(product([0,1], repeat=3))
        for i,j in combinations(range(len(unit_box_corner_vectors)), 2):
            if distance(unit_box_corner_vectors[i], unit_box_corner_vectors[j]) == 1:
                box_corner_1 = np.array(box)*unit_box_corner_vectors[i]
                box_corner_2 = np.array(box)*unit_box_corner_vectors[j]
                pts = np.array([box_corner_1, box_corner_2])
                x, y, z = pts.T
                ax.plot(x, y, z, color="r")
    if distance_cutoff_graph:
        points = np.array([xs, ys, zs]).T
        for e in distance_cutoff_graph.edges():
            point_indexes = [distance_cutoff_graph.vertex_index[e.source()], distance_cutoff_graph.vertex_index[e.target()]]
            edge_points = [points[i] for i in point_indexes]
            linestyle = "--" if cutoff_distance and distance(*edge_points) > cutoff_distance else "-"
            x, y, z = np.array(edge_points).T
            ax.plot(x, y, z, color="b", linewidth=2, linestyle=linestyle)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    plt.show()

def plot_pbc(pbc_points, box):
    xs, ys, zs = pbc_points.T
    scatter_3d(xs, ys, zs, box=box)

def plot_coms():
    models = load_models()
    for model in models.values():
        print model
        filter_model(model, ["IPR", "IPS"])
        coms = get_coms(model)
        #coms = np.array([(0,0,0)])
        xs, ys, zs = zip(*coms)
        scatter_3d(xs, ys, zs, box=model.box)

def generate_distance_cutoff_graph(points, pbc_distances, distance_cutoff, max_distance, implicit_pbc=True):
    n_pts = points.shape[0]
    g = Graph(directed=False)
    g.add_vertex(n_pts)

    g.vertex_properties["connected"] = g.new_vertex_property("bool")
    g.edge_properties["connected"] = g.new_edge_property("bool")
    for v in g.vertices():
        g.vertex_properties.connected[v] = False

    for i, j in combinations(range(n_pts), 2):
        if pbc_distances[i][j] < distance_cutoff:
            e = g.add_edge(g.vertex(i), g.vertex(j))
            g.edge_properties.connected[e] = False
    if not implicit_pbc:
        print "Filtering graph based on shortest box dimension"
        for i in range(n_pts/N_IMAGES):
            g.vertex_properties.connected[g.vertex(i)] = True
            for e in bfs_iterator(g, g.vertex(i)):
                if pbc_distances[i][g.vertex_index[e.target()]] < max_distance and g.vertex_properties.connected[e.source()]:
                    g.vertex_properties.connected[e.target()] = True
                    g.edge_properties.connected[e] = True

        g.set_filters(g.edge_properties["connected"], g.vertex_properties["connected"], )
    return g

def draw_edges_in_3D(pbc_points, distance_cutoff_graph, box, cutoff_distance=None):
    xs, ys, zs = pbc_points.T
    scatter_3d(xs, ys, zs, box=box, distance_cutoff_graph=distance_cutoff_graph, cutoff_distance=cutoff_distance)


def plot_distance_histogram(distances):
    his, bins = histogram(distances, bins = 30)
    centers = (bins[:-1]+bins[1:])/2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(centers, his)
    ax.set_xlabel("Distance (nm)")
    ax.set_ylabel("Occurrence")
    plt.show()

def plot_vertex_degree_histogram(distance_cutoff_graph):
    his, bins = vertex_hist(distance_cutoff_graph, "total")
    centers = (bins[:-1]+bins[1:])/2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(centers, his)
    ax.set_xlabel("Vertex degree")
    ax.set_ylabel("Occurrence")
    plt.show()

def main():
    cutoff_distance = 2
    models = load_models()
    average_vertex_degree = {}
    for model in models.values():
        print model
        filter_model(model, ["IPR", "IPS"])
        points = get_coms(model)#np.array([[0,0,0]])
        print "N Points: {0}".format(len(points))
        #points, pbc_point_distances = get_explicit_pbc_distances(points, model.box, )
        pbc_point_distances = get_implicit_pbc_distances(points, model.box, )

        #plot_distance_histogram(squareform(pbc_point_distances))
        #plot_coms()
        #plot_pbc(pbc_points, model.box)
        # the maximum distance to consider must be less than the smallest box dimension to avoid particles seeing themselves
        max_distance = np.min(np.array(model.box)[np.nonzero(model.box)])

        distance_cutoff_graph = generate_distance_cutoff_graph(points, pbc_point_distances, cutoff_distance, max_distance, implicit_pbc=True)
        #graph_draw(distance_cutoff_graph, vertex_text=distance_cutoff_graph.vertex_index)
        #draw_edges_in_3D(points, distance_cutoff_graph, model.box, cutoff_distance=cutoff_distance)
        #plot_vertex_degree_histogram(distance_cutoff_graph)
        average_vertex_degree[str(model)] = vertex_average(distance_cutoff_graph, "total")
    print "\n".join(["{0}: average number of neighbours within {1}nm = {2}".format(k, cutoff_distance, "{0:.2f} +/- {1:.2f}".format(*v)) for k,v in sorted(average_vertex_degree.items(), key=lambda x:int(x[0].split()[7]))])

if __name__=="__main__":
    g = Graph(directed=False)
    g.add_vertex(2)
    g.add_edge(g.vertex(0), g.vertex(1))
    for v in bfs_iterator(g, g.vertex(0)):
        print v
    #graph_draw(g, vertex_text=g.vertex_index)
    #print shortest_path(g, g.vertex(0), g.vertex(0))
    main()