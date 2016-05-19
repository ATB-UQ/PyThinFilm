import subprocess
from os.path import join, exists
from itertools import product
import numpy as np

g_density_script = '''
module load gromacs-4.6
g_density -f {gromacs_output_file} -s md.tpr -o {output_path} -n {index_file} -sl {n_divisions} -dens number << EOF
{g_density_selection}
EOF'''

make_ndx_script = '''
module load gromacs-4.6
make_ndx -f end.gro -o {index_file} << EOF
4 | 5
q
EOF
'''
SYSTEM_IDS = {"CBP":3, "IPS": 4, "IPR": 5, "BCP": 6, "IPS_IPR": 7}

END_DATA_PATH = "/mddata2/uqmstroe/MD-PhaseTransition/CBP-10nm-CBP-10nm-{0}K-slow/2200"
INIT_DATA_PATH = "/mddata2/uqmstroe/MD-PhaseTransition/CBP-10nm-CBP-10nm-{0}K-slow/2090"

TEMP_COLOURS = {"300":"k", "350":"b", "400":"g", "450":"r"}
SPECIES_LINE_STYLE = {"CBP":"-", "BCP":"--", 'IPS_IPR':"-", "BCP_20ns": "-.","BCP_50ns": "-", "CBP_20ns": ":", "CBP_50ns": "."}
N_DIVISIONS = 50
GROMACS_OUTPUT_TYPE = {"xtc": "md.xtc", "gro":"end.gro"}
DEFAULT_GROMACS_OUTPUT_TYPE = "gro"

def plot(x, y, species, temperature, figure_name, show=True, ax=None, save=False):
    import os
    if not os.environ.has_key("DISPLAY"):
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.hold(True)
    ax.plot(x, y, TEMP_COLOURS[temperature] + SPECIES_LINE_STYLE[species], label="{0} {1}K".format(species.replace("_", " | "), temperature))
    if save or show:
        ax.set_xlim(0, 20)
        plt.ylabel('Number density (nm^-3)', fontweight="bold")
        plt.xlabel('Box (nm)', fontweight="bold")
        #plt.title(figure_name[:-4].replace("_", " ").capitalize(), fontweight="bold")
        plt.legend(loc = 'upper right', prop={'size':11}, numpoints = 1, frameon = False)
    if save:
        plt.savefig(figure_name, format="eps")
    if show and os.environ.has_key("DISPLAY"):
        plt.show()
    return ax

def run_stdout_stderr_piped(script, cwd_opt=None):
    print script
    proc = subprocess.Popen(script, stderr=subprocess.PIPE, stdout=subprocess.PIPE, cwd=cwd_opt, shell=True)
    stdout, stderr = proc.communicate()
    return stdout, stderr

def run_script(script, cwd_opt=None):
    print script
    subprocess.Popen(script, cwd=cwd_opt, shell=True).wait()

def density(species, temperature, n_divisions, gromacs_output_type, data_path):
    data = data_path.format(temperature)
    output_path = join(data, "{0}_ndiv{1}_{2}.xvg".format(species, n_divisions, gromacs_output_type))
    index_file = join(data, "md.ndx")
    if not exists(index_file):
        run_script(make_ndx_script.format(index_file=index_file), cwd_opt=data)
    if not exists(output_path):
        run_script(g_density_script.format(gromacs_output_file=GROMACS_OUTPUT_TYPE[gromacs_output_type], output_path=output_path, g_density_selection=SYSTEM_IDS[species], index_file=index_file, n_divisions=n_divisions), cwd_opt=data)
    return zip(*[map(float, l.split()) for l in open(output_path).read().splitlines() if l and l[0] not in "#@"])

def plot_density(species, temperature, n_divisions, gromacs_output_type, ref_diff, show, figure_name, ax=None):
    x, ye = density(species, temperature, n_divisions, gromacs_output_type, END_DATA_PATH)
    if ref_diff:
        _, yo = get_ref_data(species, temperature, n_divisions)
        y = np.array(ye) - np.array(yo)
        diff_str = "_difference"
    else:
        y = ye
        diff_str = ""
    ax = plot(x, y, species, temperature, figure_name.format(diff_str), show=show, ax=ax, save=show)
    return ax

def plot_selected_densities(selections, n_divisions=N_DIVISIONS, gromacs_output_type=DEFAULT_GROMACS_OUTPUT_TYPE, ref_diff=False, name_prefix=""):
    ax = None
    figure_name = name_prefix + "partial_densities{0}.eps"
    for i, (temperature, species) in enumerate(selections):
        show = (i==len(selections)-1)
        ax = plot_density(species, temperature, n_divisions, gromacs_output_type, ref_diff, show, figure_name, ax=ax)

def get_ref_data(species, temperature, n_divisions):
    return density(species, temperature, n_divisions, "gro", INIT_DATA_PATH)

def density_analysis():
    selections = list(product(["300", "350", "400", "450"], ["CBP", "BCP"]))
    plot_selected_densities(selections, n_divisions=50)

    selections = list(product(["300", "350", "400", "450"], ["IPS_IPR"]))
    plot_selected_densities(selections, n_divisions=50, name_prefix="emitter_1nm_bin_width")
    plot_selected_densities(selections, n_divisions=100, name_prefix="emitter_0.5nm_bin_width_")

def integral_density_analysis():
    selections = list(product(["450"], ["CBP", "BCP"]))
    ax = None
    n_divisions = 100
    figure_name = "cumulative_densities.eps"
    perc = .95
    for i, (temperature, species) in enumerate(selections):
        _, yo = get_ref_data(species, temperature, n_divisions)
        norm_factor_ref = sum(yo)
        cumulative_sum_ref = [sum(yo[:j])/norm_factor_ref for j, _ in enumerate(yo)]
        x, ye = density(species, temperature, n_divisions, "gro", END_DATA_PATH)
        norm_factor = sum(ye)
        cumulative_sum = [sum(ye[:j])/norm_factor for j, _ in enumerate(ye)]

        show = (i==len(selections)-1)
        if species == "CBP":
            print "{0}: {1}nm".format(species, find_perc(x, cumulative_sum_ref, perc) - find_perc(x, cumulative_sum, perc))
        else:
            print "{0}: {1}nm".format(species, find_perc(x, cumulative_sum_ref, 1-perc) - find_perc(x, cumulative_sum, 1-perc))

        ax = plot(x, np.array(cumulative_sum), species + "_50ns", temperature, figure_name, show=False, ax=ax, save=False)
        ax = plot(x, np.array(cumulative_sum_ref), species + "_20ns", temperature, figure_name, show=show, ax=ax, save=show)

def find_perc(x, y, perc):
    for x, y in zip(x, y):
        if y > perc:
            return x

density_analysis()
integral_density_analysis()
