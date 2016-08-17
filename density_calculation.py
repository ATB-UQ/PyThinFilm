import os
import pmx
from plot import save_figure, add_axis_to_figure, create_figure, plot
import numpy as np
DATA_PATH = "deposition_runs/"
DEPOSITED_GRO = {
    0: os.path.join(DATA_PATH, "cbp.gro"),
    15: os.path.join(DATA_PATH, "n15_cbp.gro"),
    40: os.path.join(DATA_PATH, "n40_cbp.gro"),
    108:os.path.join(DATA_PATH, "n108_cbp.gro"),
    209:os.path.join(DATA_PATH, "n209_cbp.gro"),
    }
Z_RANGE = [3,9]
INCLUDE = ["IPR", "IPS", "CBP"]
MASS = {"C":12.011, "H":1.008, "N":14.0067, "I":192.2170}
IMAGE_FORMAT = "eps"
def load_model(model_path):
    model = pmx.Model(model_path)
    #model.nm2a()
    print "loaded: {0}".format(model)
    return model

def calculate_density(model, z_range=Z_RANGE):
    total_mass = 0
    for res in model.residues:
        if res.resname in INCLUDE:
            for atom in res.atoms:
                if z_range[0] < atom.x[2] < z_range[1]:
                    total_mass += MASS[atom.name[0]]
    volume =  model.box[0][0]*model.box[1][1]*(z_range[1] - z_range[0]) # nm^3
    return 1.6605402e-3*total_mass/volume

def single_frame_analysis_for_all():
    for n_irppy in (0, 15, 40, 108, 209):
        print n_irppy
        model = load_model(DEPOSITED_GRO[n_irppy])
        print calculate_density(model)

def plot_density(n_irppy):
    window_size = 5
    window_delta = 0.5
    min_z = 0
    max_z = 10
    density = []
    window_pos = []
    model = load_model(DEPOSITED_GRO[n_irppy])
    for left_window_pos in np.arange(min_z, max_z, window_delta):
        density.append(calculate_density(model, z_range=[left_window_pos, left_window_pos+window_size]))
        window_pos.append(left_window_pos+window_size/2.)
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot(ax, window_pos, density, color="k", xlabel="z-axis (nm)", ylabel="density (g/cm^3)", zorder=2, linewidth=1, legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./density_profile_n_irppy{0}_ws{1}_wd{2}".format(n_irppy, window_size, window_delta), image_format=IMAGE_FORMAT)
    print n_irppy
    print np.max(density)

def plot_all_density():
    for n_irppy in (0, 15, 40, 108, 209):
        plot_density(n_irppy)

if __name__=="__main__":
    #single_frame_analysis_for_all()
    #plot_density(15)
    plot_all_density()