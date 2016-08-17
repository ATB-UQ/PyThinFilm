import pmx
import os
import pickle

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
DATA_PATH = "./deposition_runs/"
DEPOSITED = {
    0: os.path.join(DATA_PATH, "cbp_only/cbp.gro"),
    15: os.path.join(DATA_PATH, "ip_15n/{0}_density.xvg"),
    40: os.path.join(DATA_PATH, "ip_40n/{0}_density.xvg"),
    108:os.path.join(DATA_PATH, "ip_108n/{0}_density.xvg"),
    209:os.path.join(DATA_PATH, "ip_209n/{0}_density.xvg"),
    }

IMAGE_FORMAT = "eps"
CONV_FACTOR = 0.001
x_lim = (0,11.5)
y_lim = (0,2)
line_styles = [(1.5, 2.5), (4,2), None]
line_width = [1.5,1,1]

TITLE_MAP = {15:"2.0 wt%",40:"5.3 wt%",108:"14.1 wt%",209:"26.3 wt%"}
def parse_xvg(xvg):
    with open(xvg) as fh:
        z, den = zip(*[map(float, l.split()) for l in fh.read().splitlines() if l and l[0] not in "#@"])
    return z, np.array(den)*CONV_FACTOR

def plot_densities(data):
    fig = create_figure((7,5))
    layout = "22{0}"
    show_x_label_in = [2,3]
    show_y_label_in = [0,2]
    for i, (n_irppy, densities) in enumerate(sorted(data.items())):
        ax = add_axis_to_figure(fig, subplot_layout=layout.format(i+1))
        for j, (species, (z, density)) in enumerate(sorted(densities.items())):
            x_label = "z (nm)" if i in show_x_label_in else None
            y_label = "density (g cm$^{-3}$)" if i in show_y_label_in else None
            plot(ax, z, density, color="k", text=TITLE_MAP[n_irppy], text_location=(0.87, 0.9), xlabel=x_label, ylabel=y_label, zorder=2, linewidth=line_width[j], xlim=x_lim, ylim=y_lim, dashes=line_styles[j])
            #plot(ax, xs, expected_distribution, color="k",zorder=1, linewidth=1, dashes=(4,2), xlim=(0,90), legend_position="upper left", legend_frame=False)
    fig.tight_layout()
    save_figure(fig, "./partial_densities", image_format=IMAGE_FORMAT)

def calculate_density_data(n_irppy):
    data = {}
    for species in ("total", "cbp", "irppy"):
        data[species] = parse_xvg(DEPOSITED[n_irppy].format(species))
    return data

def single_frame_analysis_for_all():
    density_data = {}
    for n_irppy in (15, 40, 108, 209):
        density_data[n_irppy] = calculate_density_data(n_irppy)
    plot_densities(density_data)

if __name__=="__main__":
    single_frame_analysis_for_all()
