import subprocess
from os.path import join, exists
from itertools import product
import numpy as np

g_density_script = '''
module load gromacs-4.6
g_msd -f md_nojump.xtc -s md.tpr -o {output_path} -n {index_file} << EOF
{g_diffusion_selection}
EOF'''

make_ndx_script = '''
module load gromacs-4.6
make_ndx -f end.gro -o {index_file} << EOF
4 | 5
q
EOF
'''
SYSTEM_IDS = {"CBP":3, "IPS": 4, "IPR": 5, "BCP": 6, "IPS_IPR": 7}

DATA_PATH = "/mddata2/uqmstroe/MD-PhaseTransition/CBP-10nm-CBP-10nm-{0}K-slow/2200/"

TEMP_COLOURS = {"300":"k", "350":"b", "400":"g", "450":"r"}
SPECIES_LINE_STYLE = {"CBP":"-", "BCP":"--", 'IPS_IPR':"-"}

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
        #ax.set_xlim(0, 20)
        plt.ylabel('MSD', fontweight="bold")
        plt.xlabel('x', fontweight="bold")
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

def diffusion(species, temperature, data_path):
    data = data_path.format(temperature)
    output_path = join(data, "{0}_diffusion.xvg".format(species))
    index_file = join(data, "md.ndx")
    if not exists(index_file):
        run_script(make_ndx_script.format(index_file=index_file), cwd_opt=data)
    print output_path
    if not exists(output_path):
        run_script(g_density_script.format(output_path=output_path, g_diffusion_selection=SYSTEM_IDS[species], index_file=index_file), cwd_opt=data)
    return zip(*[map(float, l.split()) for l in open(output_path).read().splitlines() if l and l[0] not in "#@"])

def plot_diffusion(species, temperature, show, figure_name, ax=None):
    x, ye = diffusion(species, temperature, DATA_PATH)
    y = ye
    ax = plot(x, y, species, temperature, figure_name, show=show, ax=ax, save=show)
    return ax

def plot_selected_diffusion(selections, name_prefix=""):
    ax = None
    figure_name = name_prefix + "diffusion{0}.eps"
    for i, (temperature, species) in enumerate(selections):
        show = (i==len(selections)-1)
        ax = plot_diffusion(species, temperature, show, figure_name, ax=ax)

def diffusion_analysis():
    selections = list(product(["300", "350", "400", "450"], ["CBP", "BCP"]))
    plot_selected_diffusion(selections)

    selections = list(product(["300", "350", "400", "450"], ["IPS_IPR"]))
    plot_selected_diffusion(selections)



diffusion_analysis()
