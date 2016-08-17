import cPickle
from os.path import exists
import numpy as np
from numpy import histogram, meshgrid
import matplotlib.pyplot as plt
from itertools import product, combinations
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, squareform
from graph_tool import Graph
from graph_tool.draw import graph_draw
from graph_tool.search import bfs_iterator
from graph_tool.stats import vertex_average, vertex_hist
from graph_tool.topology import shortest_path

import pmx
import pprint
import yaml
import traceback

MODEL_PATH_TEMPLATE = "whole_n{0}.gro"
MODEL_CACHE_TEMPLATE = "cached_modles.pickle"
USE_CACHE = False
WPC = (15, 40, 108, 209)

DIMENSIONS = 3
N_IMAGES = DIMENSIONS**2
SMALL_NUMBER = 1e-12

Z_HEIGHT_TRUNCATION = 11

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

def get_explicit_pbc_distances(points, box_lengths):
    box_x, box_y, box_z = box_lengths
    pbc_points = get_pbc_points(points, box_x, box_y, box_z)
    distances = pdist(pbc_points)
    distances = squareform(distances)
    return pbc_points, distances

def get_distance_deltas(points):
    point_compinations = list(combinations(range(points.shape[0]), 2))
    distance_deltas = np.empty((len(point_compinations), DIMENSIONS))
    for i, (j,k) in enumerate(point_compinations):
        distance_deltas[i] = points[j] - points[k]
    return distance_deltas

def get_implicit_pbc_distances(points, box):
    return get_pbc_distances_for_dims(points, box, [0,1,2])

def get_pbc_distances_for_dims(points, box, dims):
    distances = np.empty((points.shape[0], points.shape[0]))
    box_lengths = box[0][0], box[1][1], box[2][2]
    for i,j in product(range(points.shape[0]), repeat=2):
        distances[i][j] = pbc_distance_fast(points[i], points[j], box_lengths, dims=dims)
    return distances

def min_pbc_distance(point1, point2, box_lengths):
    box_x, box_y, box_z = box_lengths
    pbc_points = get_pbc_points(np.array([point1, point2]), box_x, box_y, box_z)
    return np.min(pdist(pbc_points))

def pbc_distance(point1, point2, box_lengths):
    box_x, box_y, box_z = box_lengths
    pbc_deltas = []
    for v1, v2, b in zip(point1, point2, (box_x, box_y, box_z)):
        d = v1-v2
        if abs(d) > b/2.:
            pbc_deltas.append(d+b)
        else:
            pbc_deltas.append(d)
    return np.sqrt(np.sum([v**2 for v in pbc_deltas]))

def pbc_distance_fast(point1, point2, box_lengths, dims=[0,1,2]):
    pbc_deltas = pbc_dim_deltas(point1, point2, box_lengths, dims)
    return np.sqrt(np.sum(np.power(pbc_deltas, 2)))

def pbc_dim_deltas(point1, point2, box_lengths, dims):
    pbc_deltas = np.zeros(DIMENSIONS)
    for i in dims:
        pbc_deltas[i] = point1[i]-point2[i]
        box_dim = box_lengths[i]
        if np.abs(pbc_deltas[i]) > box_dim/2.:
            pbc_deltas[i] -= box_dim*round(pbc_deltas[i]/box_dim)
    return pbc_deltas

def load_model(model_path):
    model = pmx.Model(model_path)
    #model.nm2a()
    print model
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

#TODO: finish
def scatter_3d(xs, ys, zs, box=None, distance_cutoff_graph=None, cutoff_distance=None, max_distance_vertices=None, output=None, draw_graphene=False):

    #xmin = min(xs); xmax = max(xs)
    #ymin = min(ys); ymax = max(ys)

    fig = plt.figure(figsize=(5,5.2))
    ax = fig.add_subplot(111, projection='3d')#mpl_toolkits.mplot3d.Axes3D(fig)
    ax.scatter(xs, ys, zs, alpha=0.5)#, s=300)
    if box:
        box = [box[0][0], box[1][1], box[2][2]]
        unit_box_corner_vectors = list(product([0,1], repeat=3))
        for i,j in combinations(range(len(unit_box_corner_vectors)), 2):
            if distance(unit_box_corner_vectors[i], unit_box_corner_vectors[j]) == 1:
                box_corner_1 = np.array(box)*unit_box_corner_vectors[i]
                box_corner_2 = np.array(box)*unit_box_corner_vectors[j]
                pts = np.array([box_corner_1, box_corner_2])
                x, y, z = pts.T
                zorder = 5 if (any([i > 0 for i in x]) and any([i == 0 for i in y])) else -1
                ax.plot(x, y, z, color="r", zorder=zorder)
    if distance_cutoff_graph:
        points = np.array([xs, ys, zs]).T
        for e in distance_cutoff_graph.edges():
            point_indexes = [distance_cutoff_graph.vertex_index[e.source()], distance_cutoff_graph.vertex_index[e.target()]]
            edge_points = [points[i] for i in point_indexes]
            linestyle = ":" if cutoff_distance and distance(*edge_points) > np.min(box)/2.0 else "-"
            x, y, z = np.array(edge_points).T
            ax.plot(x, y, z, color="b", linewidth=1.5, linestyle=linestyle)
    if max_distance_vertices:
        for vertex_pair in max_distance_vertices:
            point_indexes = [distance_cutoff_graph.vertex_index[vertex_pair[0]], distance_cutoff_graph.vertex_index[vertex_pair[1]]]
            edge_points = [points[i] for i in point_indexes]
            linestyle = "-"#"--" if cutoff_distance and distance(*edge_points) > np.min(box)/2.0 else "-"
            x, y, z = np.array(edge_points).T
            ax.plot(x, y, z, color="g", linewidth=3, linestyle=linestyle)
    if draw_graphene:
        nx, ny = (20, 20)
        x = np.linspace(0, box[0], nx)
        y = np.linspace(0, box[1], ny)
        X, Y = np.meshgrid(x, y)
        Z = np.zeros(X.shape)
        ax.plot_wireframe(X, Y, Z, zorder=-1, alpha=0.5, color="#BEC8CC")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    #print "\n".join(n for n in dir(ax) if "tick" in n)
    ticks = range(0,10,2)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    padding = 0.2
    #max_val = np.max(np.concatenate([xs, ys, zs]))
    ax.set_xlim((-padding, box[0] + padding ))
    ax.set_ylim((-padding, box[1] + padding))
    ax.set_zlim((-padding, box[2] + padding))
    fig.subplots_adjust(top=1, left=0, right=1, bottom=0)
    ax.dist = 10.8
    if output:
        fig.savefig(output)
    else:
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

def generate_distance_cutoff_graph(points, pbc_distances, distance_cutoff, box_lengths=None, implicit_pbc=True):
    # if box_lengths is not None, check whether point distances are greater than box lengths
    box_length_exeeded = False
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
        #print "Filtering graph based on box dimensions"
        for i in range(n_pts/N_IMAGES):
            g.vertex_properties.connected[g.vertex(i)] = True
            for e in bfs_iterator(g, g.vertex(i)):
                if g.vertex_properties.connected[e.source()]:
                    too_long = False
                    if box_lengths is not None:
                        distance_deltas = points[i] - points[g.vertex_index[e.target()]]
                        for bl, dim_delta in zip(box_lengths, distance_deltas):
                            too_long = dim_delta > bl
                            if too_long:
                                break
                        if not box_length_exeeded and too_long:
                            box_length_exeeded = True
                    if not too_long:
                        g.vertex_properties.connected[e.target()] = True
                        g.edge_properties.connected[e] = True

        g.set_filters(g.edge_properties["connected"], g.vertex_properties["connected"], )
    if box_lengths:
        return g, box_length_exeeded
    return g

def draw_edges_in_3D(pbc_points, distance_cutoff_graph, box, cutoff_distance=None, max_distance_vertices=None, output=None):
    xs, ys, zs = pbc_points.T
    scatter_3d(xs, ys, zs, box=box, distance_cutoff_graph=distance_cutoff_graph, cutoff_distance=cutoff_distance, max_distance_vertices=max_distance_vertices, output=output)

def pbc_path_distances(g, source, target, points, box_lengths, dims):
    '''Calculate the through space distance in the dimensions provided by dim (x=0, y=1, z=2). Due to periodic images,
    distance must be calculated from deltas between the nodes along a given path and not simply from the end-points.'''
    if source == target:
        return 0.0
    _, path_edges = shortest_path(g, source, target)
    path_deltas = np.empty((len(path_edges), 3))
    for i, e in enumerate(path_edges):
        path_deltas[i] = pbc_dim_deltas(points[g.vertex_index[e.source()]], points[g.vertex_index[e.target()]], box_lengths, dims)
    # calculate cumulative path distances in x,y,z;
    return total_distance_from_deltas(path_deltas)

def total_distance_from_deltas(path_deltas):
    distance_vector = np.sum(path_deltas, axis=0)
    return np.sqrt(np.sum(np.power(distance_vector, 2)))

def calc_path_distances_from_path_indexes(path_indexes, points, box_lengths):
    n_edges = len(path_indexes)-1
    path_deltas = np.empty((n_edges, 3))
    for i, vertex_index_i in enumerate(path_indexes):
        j = i+1
        if j > n_edges:
            break
        path_deltas[i] = pbc_dim_deltas(points[vertex_index_i], points[path_indexes[j]], box_lengths, [0,1,2])
    # calculate cumulative path distances in x,y,z;
    return total_distance_from_deltas(path_deltas)

def is_pbc_circular(g, node_cluster, points, box_lengths, distance_cutoff, plot=False):
    cluster_points = np.array([points[g.vertex_index[v]] for v in node_cluster])
    pbc_cluster_points, pbc_point_distances = get_explicit_pbc_distances(cluster_points, box_lengths)
    pbc_cluster_graph, box_length_exceeded = generate_distance_cutoff_graph(pbc_cluster_points, pbc_point_distances, distance_cutoff, box_lengths=box_lengths, implicit_pbc=False)
    if plot:
        box = [[box_lengths[0],0,0],[0,box_lengths[1],0],[0,0,box_lengths[2]]]
        draw_edges_in_3D(pbc_cluster_points, pbc_cluster_graph, box=box, cutoff_distance=distance_cutoff)

    return box_length_exceeded

def calc_max_connected_distances(g, distances, points, box, distance_cutoff, dim=[0, 1, 2], plot=False, return_vertices=False, plot_connected_clusters=False):
    box_lengths = box[0][0], box[1][1], box[2][2]
    # generate list of connected node clusters
    counted = [] # keep track of which nodes have already been considered
    connected_node_clusters = []
    for v in g.vertices():
        if v in counted:
            continue
        counted.append(v) # keep track of which nodes have already been considered
        connected_nodes = [v]
        for e in bfs_iterator(g, v):
            connected_nodes.append(e.target())
            counted.append(e.target()) # keep track of which nodes have already been considered
        connected_node_clusters.append(connected_nodes)
    # for each cluster of connected_nodes, calculate the through space pair-wise distances
    max_through_space_distance = 0
    max_distance_vertices = None
    for node_cluster in connected_node_clusters:
        # check if node_cluster is circular in pbc, if so, set all distances to inf
        # can ignore case where only dimension 2 (z) is being considered
        if dim != [2] and len(node_cluster) > 1 and is_pbc_circular(g, node_cluster, points, box_lengths, distance_cutoff, plot=plot_connected_clusters):
            print "Circular cluster found, all distances will be set to inf"
            max_through_space_distance = np.Inf
            max_distance_vertices = (node_cluster[0], node_cluster[1])
            break
        for i,j in combinations(range(len(node_cluster)), 2):
            path_length = pbc_path_distances(g, node_cluster[i], node_cluster[j], points, box_lengths, dim)# distances[g.vertex_index[node_cluster[i]]][g.vertex_index[node_cluster[j]]]
            if path_length > max_through_space_distance:
                max_through_space_distance = path_length
                max_distance_vertices = (node_cluster[i], node_cluster[j])
    stats = {"max_through_space_distance": pretty_string(max_through_space_distance)}
                #average_through_space_distance=pretty_string(np.mean(all_connected_distances)),
                #std_through_space_distance=pretty_string(np.std(all_connected_distances)),
                #median_through_space_distance=pretty_string(np.median(all_connected_distances))
                #)
    if return_vertices:
        return stats, max_distance_vertices
    else:
        return stats

def plot_distance_histogram(distances):
    his, bins = histogram(distances, bins = 100)
    centers = (bins[:-1]+bins[1:])/2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(centers, his)
    ax.set_xlabel("Distance (nm)")
    ax.set_ylabel("Occurrence")
    plt.show()

def vertex_degree_histogram(distance_cutoff_graph, show=True):
    his, bins = vertex_hist(distance_cutoff_graph, "total")
    centers = (bins[:-1]+bins[1:])/2
    if show:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(centers, his)
        ax.set_xlabel("Vertex degree")
        ax.set_ylabel("Occurrence")
        plt.show()
    return his[0]

def pretty_string(x):
    return "{0:.2f}".format(float(x))


def get_edges_as_point_indexes(distance_cutoff_graph, exclude_pbc_image_crossing=False, points=None, cutoff_distance=None):
    edges = []
    for e in distance_cutoff_graph.edges():
        es = distance_cutoff_graph.vertex_index[e.source()]
        et = distance_cutoff_graph.vertex_index[e.target()]
        if exclude_pbc_image_crossing and distance(points[es], points[et]) > cutoff_distance:
            continue
        edges.append([es, et])

    return edges

def main():
    from generate_pymol_images import draw_connectivities, init, close
    init()

    cutoff_distances = [1.25, 1.5, 1.75, 2.0]
    models = load_models()
    stats = {}
    for cutoff_distance in cutoff_distances:
        for model in models.values():
            print model
            filter_model(model, ["IPR", "IPS"])
            points = get_coms(model)#np.array([[0,0,0]])
            model_id = "{0}nm_cutoff_{1}".format(cutoff_distance, len(points))
            stats[model_id] = {}
            print "N Points: {0}".format(len(points))
            stats[model_id]["N"] = len(points)

            pbc_point_distances = get_implicit_pbc_distances(points, model.box, )

            #plot_distance_histogram(squareform(pbc_point_distances))
            #plot_coms()
            #plot_pbc(pbc_points, model.box)

            distance_cutoff_graph = generate_distance_cutoff_graph(points, pbc_point_distances, cutoff_distance, implicit_pbc=True)
            try:
                graph_draw(distance_cutoff_graph,
                           #vertex_text=distance_cutoff_graph.vertex_index,
                           output="images/{0}_connectivity_graph.pdf".format(model_id),
                           vertex_fill_color="#1938FC",
                           vertex_pen_width=0.2,
                           output_size=(300,300),
                           edge_pen_width=1.5,
                           edge_color="#000000",
                           vertex_size=8,
                           )
            except:
                print traceback.format_exc()
            box = model.box
            box[2][2] = Z_HEIGHT_TRUNCATION
            draw_edges_in_3D(points, distance_cutoff_graph, box, cutoff_distance=cutoff_distance, output="images/{0}_3d_connectivities.pdf".format(model_id))
            edges = get_edges_as_point_indexes(distance_cutoff_graph, exclude_pbc_image_crossing=True, points=points, cutoff_distance=cutoff_distance)
            draw_connectivities("whole_n{0}".format(len(points)), points, edges, cutoff_distance, box_lengths=10*np.array([box[0][0], box[1][1], box[2][2]]))

            #N_without_edge = vertex_degree_histogram(distance_cutoff_graph, show=False)
            #mu_deg, sigma_deg = vertex_average(distance_cutoff_graph, "total")
            #stats[model_id]["average_vertex_degree"] = pretty_string(mu_deg)
            #stats[model_id]["std_vertex_degree"] = pretty_string(sigma_deg)
            #stats[model_id]["n_without_neighbour"] = pretty_string(N_without_edge)
            #stats[model_id]["perc_without_neighbour"] = pretty_string(100.*N_without_edge/float(stats[model_id]["N"]))
            #connected_distance_stats, max_distance_vertices = calc_max_connected_distances(distance_cutoff_graph, pbc_point_distances, points, model.box, cutoff_distance, return_vertices=True)#, plot_connected_clusters=True)
            #xy_connected_distance_stats, xy_max_distance_vertices = calc_max_connected_distances(distance_cutoff_graph, pbc_point_distances, points, model.box, cutoff_distance, dim=[0, 1], return_vertices=True)
            #z_connected_distance_stats, z_max_distance_vertices = calc_max_connected_distances(distance_cutoff_graph, pbc_point_distances, points, model.box, cutoff_distance, dim=[2], return_vertices=True)
            #draw_edges_in_3D(points, distance_cutoff_graph, model.box, cutoff_distance=cutoff_distance, max_distance_vertices=[max_distance_vertices, z_max_distance_vertices, xy_max_distance_vertices])
            #stats[model_id].update({
            #    "3d": connected_distance_stats,
            #    "vertical": z_connected_distance_stats,
            #    "lateral": xy_connected_distance_stats,
            #    })
    pprint.pprint(stats)
    with open("path_stats.yml", "w") as fh:
        yaml.dump(stats, fh)
    close()

if __name__=="__main__":
    g = Graph(directed=False)
    g.add_vertex(4)
    g.add_edge(g.vertex(0), g.vertex(1))
    g.add_edge(g.vertex(1), g.vertex(2))
    g.add_edge(g.vertex(2), g.vertex(3))
    g.add_edge(g.vertex(3), g.vertex(0))
    #for path in all_paths(g, 0, 0):
    #    print path
    #graph_draw(g, vertex_text=g.vertex_index)
    #print shortest_path(g, g.vertex(0), g.vertex(0))

    main()