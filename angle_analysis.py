import pmx
import os
import pickle
import time
import pylab
os.environ["GMX_DLL"]="/home/uqbcaron/PROGRAMMING_PROJECTS/CPP/gromacs-4.0.7/build/lib/"
from pmx.xtc import Trajectory
import numpy as np
from plot import create_figure, add_axis_to_figure, plot, save_figure
#DEPOSITED_GRO = {
#    15: "/mddata2/uqmstroe/Claires_deposition_data/2wpc/999/end.gro",
#    40: "/mddata2/uqmstroe/Claires_deposition_data/6wpc/300K5ns/end5ns300K.gro",
#    108:"/mddata2/uqmstroe/Claires_deposition_data/midRatio15wpc/eq300K/300K5ns",
#    209:"/mddata2/uqmstroe/Claires_deposition_data/highRatio/300K5ns/end5ns300K.gro",
#    }
#DEPOSITED_TRAJ = {
#    15: "/mddata2/uqmstroe/Claires_deposition_data/2wpc/999/md.xtc",
#    40: "/mddata2/uqmstroe/Claires_deposition_data/6wpc/300K5ns/md.xtc",
#    108: "/mddata2/uqmstroe/Claires_deposition_data/midRatio15wpc/eq300K/300K5ns/md.xtc",
#    209: "/mddata2/uqmstroe/Claires_deposition_data/highRatio/300K5ns/md.xtc",
#    }
DATA_PATH = "deposition_runs/"
DEPOSITED_GRO = {
    0: os.path.join(DATA_PATH, "cbp.gro"),
    15: os.path.join(DATA_PATH, "n15_cbp.gro"),
    40: os.path.join(DATA_PATH, "n40_cbp.gro"),
    108:os.path.join(DATA_PATH, "n108_cbp.gro"),
    209:os.path.join(DATA_PATH, "n209_cbp.gro"),
    }
DEPOSITED_TRAJ = {
    0: os.path.join(DATA_PATH, "cbp_only/cbp.xtc"),
    15: os.path.join(DATA_PATH, "ip_15n/whole_n15.xtc"),
    40: os.path.join(DATA_PATH, "ip_40n/whole_n40.xtc"),
    108:os.path.join(DATA_PATH, "ip_108n/whole_n108.xtc"),
    209:os.path.join(DATA_PATH, "ip_209n/whole_n209.xtc"),
    }

FRAME_COUNT = {
    0: 100,
    15: 200,
    40: 10000,
    108:13000,
    209:10000,
    }

N_FRAMES = 200
IMAGE_FORMAT = "eps"
N_BINS = 100
REFERENCE_AXIS = np.array([0.,0.,1.])
CACHE_TEMPLATE = "{n_irppy}_{keep_every_ith}.pickle"

def load_model(model_path):
    model = pmx.Model(model_path)
    #model.nm2a()
    print "loaded: {0}".format(model)
    return model

def irppy_c3_axis(atoms):
    N1, N2, N3 = [a for a in atoms if "N" in a.name]
    #print N1.name, N2.name, N3.name
    N1_coords, N2_coords, N3_coords = np.array(N1.x), np.array(N2.x), np.array(N3.x)
    c3_axis = np.cross(N3_coords - N1_coords, N2_coords - N1_coords)
    return c3_axis

def cbp_long_axis(atoms, use_nitrogens=True):
    if use_nitrogens:
        atom1, atom2 = [a for a in atoms if "N" in a.name]
    else:
        atom1, atom2 = [a for a in atoms if a.name == "C1"][0], [a for a in atoms if a.name == "C7"][0]
    atom1_coords, atom2_coords = np.array(atom1.x), np.array(atom2.x)
    molecular_axis = (atom1_coords - atom2_coords)

    return molecular_axis/np.linalg.norm(molecular_axis)

# Compute the angle between two vectors (numpy.array), in radians
def angle(vector1, vector2):
    return np.arccos( np.dot(vector1, vector2) / ( np.linalg.norm(vector1) * np.linalg.norm(vector2) )  )

def angle_degrees(vector1, vector2, zero_to_90=False):
    if not zero_to_90:
        transform_data = lambda x:x
    else:
        transform_data = angles_between_0_and_90
    return transform_data(angle(vector1, vector2) * 180. / np.pi)

def angles_between_0_and_90(x):
    return -np.abs((x-90)) + 90

def uniform_random_distribution(zero_to_90, xs=None):
    end_point = np.pi / 2. if zero_to_90 else np.pi
    amplitude = 1.0 if zero_to_90 else 0.5
    xs = np.linspace(0, end_point, 1000) if xs is None else xs
    expected_distribution = amplitude * np.sin(xs) / (180 / np.pi)
    return np.rad2deg(xs), expected_distribution

def plot_hist(values, n_bins=N_BINS, label="Ir(ppy)3", xlabel="C$_3$ axis (deg)", zero_to_90=False, external_ax=None):
    xs, expected_distribution = uniform_random_distribution(zero_to_90)
    his, bins = np.histogram(values, bins = n_bins, normed=True)
    centers = (bins[:-1]+bins[1:])/2
    xlim = (0,90) if zero_to_90 else (0,180)
    if external_ax is None:
        fig = create_figure((4,3))
        ax = add_axis_to_figure(fig)
    else:
        ax=external_ax
    plot(ax, centers, his, color="b", xlabel=xlabel, ylabel="probability", zorder=2, linewidth=1.5)
    plot(ax, xs, expected_distribution, color="k",zorder=1, linewidth=1.5, dashes=(4,2), xlim=xlim, legend_position="upper left", legend_frame=False)
    if external_ax is None:
        fig.tight_layout()
        save_figure(fig, "./{0}_angle_distribution".format(label), image_format=IMAGE_FORMAT)

def load_frame_data(n_irppy, molecular_axis_definition=irppy_c3_axis, molecules_to_include=None, zero_to_90=False, com_z_include=None, use_cache=True):
    keep_every_ith = max([FRAME_COUNT[n_irppy]/N_FRAMES, 1])
    cache_file = CACHE_TEMPLATE.format(n_irppy=n_irppy, keep_every_ith=keep_every_ith)
    if use_cache and os.path.exists(cache_file):
        with open(cache_file) as fh:
            return pickle.load(fh)
    model = load_model(DEPOSITED_GRO[n_irppy])
    traj = Trajectory(DEPOSITED_TRAJ[n_irppy])
    molecular_axis_vs_reference_angles = []
    for i, _ in enumerate(traj):
        traj.update(model)
        if i % keep_every_ith != 0:
            continue
        print i
        for molecule in model.residues:
            if molecule.resname in molecules_to_include and com_z_include[0] < molecule.com(vector_only=True)[2] < com_z_include[1]:
                molecular_axis = molecular_axis_definition(molecule.atoms)
                molecular_axis_vs_reference_angles.append( angle_degrees(molecular_axis, REFERENCE_AXIS, zero_to_90=zero_to_90) )
        if not molecular_axis_vs_reference_angles:
            return []
    if use_cache:
        with open(cache_file, "w") as fh:
            pickle.dump(molecular_axis_vs_reference_angles, fh)
    return molecular_axis_vs_reference_angles

def trajectory_analysis(n_irppy, plot=False, molecular_axis_definition=None, molecules_to_include=None, zero_to_90=False, com_z_include=[0,100], use_cache=True):
    molecular_axis_vs_reference_angles = load_frame_data(n_irppy,
        molecular_axis_definition=molecular_axis_definition,
        molecules_to_include=molecules_to_include, com_z_include=com_z_include, zero_to_90=zero_to_90, use_cache=use_cache)
    if plot:
        plot_hist(molecular_axis_vs_reference_angles)
    return molecular_axis_vs_reference_angles

def single_frame_analysis(n_irppy, plot=False, molecules_to_include=None, molecular_axis_definition=None, com_z_include=[0, 1000], inverse_selection=False, zero_to_90=False):
    if inverse_selection:
        select = lambda x:not x
    else:
        select = lambda x:x

    model = load_model(DEPOSITED_GRO[n_irppy])
    print "Selecting molecules"
    molecules = []
    for molecule in model.residues:
        if molecule.resname in molecules_to_include and select(com_z_include[0] < molecule.com(vector_only=True)[2] < com_z_include[1]):
            molecules.append(molecule)
    print "{0} molecules found".format(len(molecules))
    molecular_axis_vs_reference_angles = []
    for mol in molecules:
        #print irppy.id
        molecular_axis = molecular_axis_definition(mol.atoms)
        molecular_axis_vs_reference_angles.append( angle_degrees(molecular_axis, REFERENCE_AXIS, zero_to_90=zero_to_90) )
    if plot:
        plot_hist(molecular_axis_vs_reference_angles)
    return molecular_axis_vs_reference_angles

def single_frame_analysis_for_all(species, ax=None):
    angles = []
    for n_irppy in (15, 40, 108, 209):
        angles.extend(single_frame_analysis(n_irppy), molecules_to_include=species)
    plot_hist(angles, ax=ax)

def irppy_traj_analysis_for_all(ax=None):
    angles = []
    for n_irppy in (15, 40, 108, 209):
        angles.extend(trajectory_analysis(n_irppy, molecular_axis_definition=irppy_c3_axis, molecules_to_include=["IPR", "IPS"], zero_to_90=False))
    plot_hist(angles, external_ax=ax)

def calc_difference_area(xs, expected, observed):
    absolute_difference = abs(np.array(observed) - np.array(expected))
    return np.trapz(absolute_difference, xs)

def cbp_single_frame_analysis():
    angles = single_frame_analysis(0, molecules_to_include=["CBP"], molecular_axis_definition=cbp_long_axis, com_z_include=[4, 8], inverse_selection=False, zero_to_90=True)
    plot_hist(angles, label="CBP", xlabel="CBP long axis (deg)", zero_to_90=True)

def cbp_traj_analysis():
    angles = trajectory_analysis(0, molecular_axis_definition=cbp_long_axis, molecules_to_include=["CBP"], zero_to_90=True, com_z_include=[4,6])
    plot_hist(angles, label="CBP", xlabel="CBP long axis (deg)", zero_to_90=True)

def cbp_orientation_along_z_traj():
    window_size = 2
    window_delta = 1
    min_z = 0
    max_z = 11
    mean_angles_along_z = []
    window_centers = []
    for left_window_pos in np.arange(min_z, max_z, window_delta):
        molecular_axis_vs_reference_angles = trajectory_analysis(0,
            molecular_axis_definition=cbp_long_axis,
            molecules_to_include=["CBP"],
            zero_to_90=True,
            com_z_include=[left_window_pos, left_window_pos+window_size],
            use_cache=False)
        if not molecular_axis_vs_reference_angles:
            continue
        mean_angles_along_z.append(np.mean(molecular_axis_vs_reference_angles))
        window_centers.append(left_window_pos+window_size/2.)
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot(ax, window_centers, np.abs(mean_angles_along_z-np.rad2deg(1)), color="k", xlabel="z-axis (nm)", ylabel="orientational order (deg)", zorder=2, linewidth=1, legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./mean_angle_vs_z_distance_ws{0}_wd{1}".format(window_size, window_delta), image_format=IMAGE_FORMAT)

def cbp_orientation_along_z(n_irppy):
    window_size = 0.5
    window_delta = 0.1
    min_z = 0
    max_z = 10.5
    model = load_model(DEPOSITED_GRO[n_irppy])
    print "Selecting molecules"
    all_cbp_molecules = []
    for molecule in model.residues:
        if molecule.resname == "CBP":
            all_cbp_molecules.append(molecule)
    print "{0} molecules found".format(len(all_cbp_molecules))
    mean_angles_along_z = []
    window_pos = []
    for left_window_pos in np.arange(min_z, max_z, window_delta):
        molecules = [molecule for molecule in all_cbp_molecules if left_window_pos < molecule.com(vector_only=True)[2] < left_window_pos+window_size]
        molecular_axis_vs_reference_angles = []
        print left_window_pos
        if not molecules:
            continue
        for mol in molecules:
            molecular_axis = cbp_long_axis(mol.atoms)
            molecular_axis_vs_reference_angles.append( angle_degrees(molecular_axis, REFERENCE_AXIS, zero_to_90=True) )
        mean_angles_along_z.append(np.mean(molecular_axis_vs_reference_angles))
        window_pos.append(left_window_pos)
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot(ax, window_pos, (mean_angles_along_z-np.rad2deg(1))/(90-np.rad2deg(1)), color="k", xlabel="z-axis (nm)", ylim=None, xlim=(min_z, max_z+window_size/2.), ylabel="orientational order (deg)", zorder=2, linewidth=1, legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./mean_angle_vs_z_distance_n_irppy{0}_ws{1}_wd{2}".format(n_irppy, window_size, window_delta), image_format=IMAGE_FORMAT)
    return window_pos, (mean_angles_along_z-np.rad2deg(1))/(90-np.rad2deg(1)), min_z, max_z+window_size/2.

def mean_cbp_orientation_along_z(n_irppy):
    window_size = 0.5
    window_delta = 0.1
    min_z = 0
    max_z = 10.5
    model = load_model(DEPOSITED_GRO[n_irppy])
    print "Selecting molecules"
    all_cbp_molecules = []
    for molecule in model.residues:
        if molecule.resname == "CBP":
            all_cbp_molecules.append(molecule)
    print "{0} molecules found".format(len(all_cbp_molecules))
    mean_angles_along_z = []
    window_pos = []
    for left_window_pos in np.arange(min_z, max_z, window_delta):
        molecules = [molecule for molecule in all_cbp_molecules if left_window_pos < molecule.com(vector_only=True)[2] < left_window_pos+window_size]
        molecular_axis_vs_reference_angles = []
        print left_window_pos
        if not molecules:
            continue
        for mol in molecules:
            molecular_axis = cbp_long_axis(mol.atoms)
            molecular_axis_vs_reference_angles.append( angle_degrees(molecular_axis, REFERENCE_AXIS, zero_to_90=True) )
        mean_angles_along_z.append(np.mean(molecular_axis_vs_reference_angles))
        window_pos.append(left_window_pos)
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot(ax, window_pos, mean_angles_along_z, color="k", xlabel="z-axis (nm)", ylim=None, xlim=(min_z, max_z+window_size/2.), ylabel="orientational order (deg)", zorder=2, linewidth=1, legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./mean_angle_vs_z_distance_n_irppy{0}_ws{1}_wd{2}".format(n_irppy, window_size, window_delta), image_format=IMAGE_FORMAT)
    return window_pos, mean_angles_along_z, min_z, max_z+window_size/2.


def cbp_orientation_for_all(ax_external=None):
    if ax_external is None:
        fig = create_figure((4,3))
        ax = add_axis_to_figure(fig)
    else:
        ax = ax_external
    dashes = [None, (1.5,2.5)]#(3.5,1.5)]
    line_widths = [1.5, 1.75]
    line_styles = ["-", "-"]
    color = ["b", "k"]
    for i, n_irppy in enumerate([0, 40]):#, 15, 40, 108, 209):
        window_pos, mean_angles_along_z, min_z, max_z = cbp_orientation_along_z(n_irppy)
        plot(ax, window_pos, mean_angles_along_z, color=color[i], xlabel="z-axis (nm)", ylim=None, xlim=(min_z, max_z), ylabel="orientation order", zorder=i, linewidth=line_widths[i], legend_frame=False, line_style=line_styles[i], dashes=dashes[i])
    #plot(ax, [window_pos[0], window_pos[-1]], [np.rad2deg(1), np.rad2deg(1)], zorder=2, linewidth=1, legend_frame=False, dashes=(3.5,1.5))
    if ax_external is None:
        fig.tight_layout()
        save_figure(fig, "./mean_angle_vs_z_distance_n_irppy", image_format=IMAGE_FORMAT)
if __name__=="__main__":
    #fig = create_figure((4,6))
    #irppy_single_frame_analysis_for_all()
    #irppy_traj_analysis_for_all(add_axis_to_figure(fig, subplot_layout=212))
    #irppy_traj_analysis_for_all()
    #cbp_single_frame_analysis()
    #cbp_traj_analysis()
    #cbp_orientation_for_all(add_axis_to_figure(fig, subplot_layout=211))
    cbp_orientation_for_all()
    #fig.tight_layout()
    #save_figure(fig, "./angle_distribution_analysis", image_format=IMAGE_FORMAT)
