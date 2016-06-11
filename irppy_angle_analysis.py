import pmx
import os
import pickle
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
DATA_PATH = "/mddata2/uqmstroe/OLED_data/deposited_irppy/"
DEPOSITED_GRO = {
    15: os.path.join(DATA_PATH, "ip_15n/whole_n15.gro"),
    40: os.path.join(DATA_PATH, "ip_40n/whole_n40.gro"),
    108:os.path.join(DATA_PATH, "ip_108n/whole_n108.gro"),
    209:os.path.join(DATA_PATH, "ip_209n/whole_n209.gro"),
    }
DEPOSITED_TRAJ = {
    15: os.path.join(DATA_PATH, "ip_15n/whole_n15.xtc"),
    40: os.path.join(DATA_PATH, "ip_40n/whole_n40.xtc"),
    108:os.path.join(DATA_PATH, "ip_108n/whole_n108.xtc"),
    209:os.path.join(DATA_PATH, "ip_209n/whole_n209.xtc"),
    }

FRAME_COUNT = {
    15: 200,
    40: 10000,
    108:13000,
    209:10000,
    }

N_FRAMES = 500

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

# Compute the angle between two vectors (numpy.array), in radians
def angle(vector1, vector2):
    return np.arccos( np.dot(vector1, vector2) / ( np.linalg.norm(vector1) * np.linalg.norm(vector2) )  )

def angle_degrees(vector1, vector2):
    return angle(vector1, vector2) * 180. / np.pi

def plot_hist(values, n_bins=100):
    n = len(values)
    xs = np.linspace(0, np.pi, 100)
    expected_distribution = 0.5*np.sin(xs)/(180/np.pi)
    xs *= (180/np.pi)
    his, bins = np.histogram(values, bins = n_bins, normed=True)
    centers = (bins[:-1]+bins[1:])/2
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot(ax, centers, his, color="k", label="IrPPY3", xlabel="C3 axis angle (deg)", ylabel="Probability", zorder=2, linewidth=1)
    plot(ax, xs, expected_distribution, color="k", label="Random", zorder=1, linewidth=1, line_style="--", xlim=(0,180), legend_position="upper left", legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./c3_angle_distribution", image_format="pdf")
    fig.show()

def load_frame_data(n_irppy):
    keep_every_ith = max([FRAME_COUNT[n_irppy]/N_FRAMES, 1])
    cache_file = CACHE_TEMPLATE.format(n_irppy=n_irppy, keep_every_ith=keep_every_ith)
    if os.path.exists(cache_file):
        with open(cache_file) as fh:
            return pickle.load(fh)
    model = load_model(DEPOSITED_GRO[n_irppy])
    traj = Trajectory(DEPOSITED_TRAJ[n_irppy])
    c3_axis_vs_z_angles = []
    for i, _ in enumerate(traj):
        traj.update(model)
        if i % keep_every_ith != 0:
            continue
        print i
        for irppy in model.residues:
            #print irppy.id
            c3_axis = irppy_c3_axis(irppy.atoms)
            z_axis = np.array([0.,0.,1.])
            c3_axis_vs_z_angles.append( angle_degrees(c3_axis, z_axis) )
    with open(cache_file, "w") as fh:
        pickle.dump(c3_axis_vs_z_angles, fh)
    return c3_axis_vs_z_angles

def trajectory_analysis(n_irppy, plot=False):
    c3_axis_vs_z_angles = load_frame_data(n_irppy)
    if plot:
        plot_hist(c3_axis_vs_z_angles)
    return c3_axis_vs_z_angles

def single_frame_analysis(n_irppy, plot=False):
    model = load_model(DEPOSITED_GRO[n_irppy])
    print "Selecting irppy molecules"
    irppy_molecules = []
    for molecule in model.residues:
        if molecule.resname in ["IPR", "IPS"] and molecule.com(vector_only=True)[2] > 0:
            irppy_molecules.append(molecule)
    print "{0} irppy molecules found".format(len(irppy_molecules))
    c3_axis_vs_z_angles = []
    for irppy in irppy_molecules:
        #print irppy.id
        c3_axis = irppy_c3_axis(irppy.atoms)
        z_axis = np.array([0.,0.,1.])
        c3_axis_vs_z_angles.append( angle_degrees(c3_axis, z_axis) )
    if plot:
        plot_hist(c3_axis_vs_z_angles)
    return c3_axis_vs_z_angles

def single_frame_analysis_for_all():
    angles = []
    for n_irppy in (15, 40, 108, 209):
        angles.extend(single_frame_analysis(n_irppy))
    plot_hist(angles)

def traj_analysis_for_all():
    angles = []
    for n_irppy in (15, 40, 108, 209):
        angles.extend(trajectory_analysis(n_irppy))
    plot_hist(angles)
if __name__=="__main__":
    #single_frame_analysis_for_all()
    traj_analysis_for_all()
