from plot import create_figure, add_axis_to_figure, plot, save_figure
from matplotlib import pyplot
DEPOSITED_DATA_PATH = "/home/uqmstroe/workspace/RDF_Fun/n{0}/deposited_data/rdf_IP_IP.xvg"
RANDOM_DATA_PATH = "/home/uqmstroe/workspace/RDF_Fun/n{0}/rdf_IPD_IPD.xvg"

def plot_single_rdf(deposited, random, name):
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot_to_axis(ax, deposited, random)
    fig.tight_layout()
    save_figure(fig, "./rdf_{0}".format(name), image_format="pdf")

def plot_to_axis(ax, deposited, random):
    plot(ax, deposited[0], deposited[1], color="k", label="IrPPY3", xlabel="r (nm)", ylabel="g(r)", zorder=2, linewidth=1)
    plot(ax,  random[0], random[1], color="k", label="Random", zorder=1, linewidth=1, xlim=(0,4), legend_position="upper right", dashes=(4,2))

def plot_all_on_single_fig():
    fig = create_figure((12,3))
    ax_list = [add_axis_to_figure(fig, subplot_layout=int("14{0}".format(i+1))) for i in range(4)]
    for ax, n_irppy in zip(ax_list, (15, 40, 108, 209)):
        deposited = parse_xvg(DEPOSITED_DATA_PATH.format(n_irppy))
        random = parse_xvg(RANDOM_DATA_PATH.format(n_irppy))
        plot_to_axis(ax, deposited, random)
    fig.tight_layout()
    save_figure(fig, "./rdf_combined", image_format="pdf")
    fig.show()

def parse_xvg(data_path):
    with open(data_path) as fh:
        xvg_str = fh.read()
    xs, ys = zip(*[map(float, l.split()) for l in xvg_str.splitlines() if l[0] not in ["#", "@"]])
    return xs, ys

def plot_all():
    for n_irppy in (15, 40, 108, 209):
        deposited = parse_xvg(DEPOSITED_DATA_PATH.format(n_irppy))
        random = parse_xvg(RANDOM_DATA_PATH.format(n_irppy))
        plot_single_rdf(deposited, random, "irppy_n{0}".format(n_irppy))

if __name__=="__main__":
    plot_all_on_single_fig()
    plot_all()
