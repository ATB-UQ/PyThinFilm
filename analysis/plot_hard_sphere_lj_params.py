import sys
from itertools import product
import colorsys
from os.path import basename
import os
if not os.environ.has_key("DISPLAY"):
    import matplotlib
    matplotlib.use("Agg")
import pylab
import numpy as np

def create_figure(figsize=(10,8)):
    fig = pylab.figure(figsize=figsize)
    fig.hold = True
    return fig

def add_axis_to_figure(fig, subplot_layout=111, sharex=None, sharey=None):
    return fig.add_subplot(subplot_layout, sharex=sharex, sharey=sharey)

def save_figure(fig_name):
    file_name = '{0}.png'.format(fig_name)
    print "saving to file: {0}".format(file_name)
    pylab.savefig(file_name)

def plot(ax, xs, ys, xlabel=None, ylabel=None, label=None, linewidth=1, alpha=1):
    marker_size = 10
    width = 1.25 # table borders and ticks
    tick_width = 0.75
    tick_length = 4
    font_size = 12
    font_weight = "bold"

    symbols = [u'o', u'v', u'^', u'<', u'>', u's', u'p', u'h', u'H', u'D', u'd']
    symbol = None#symbols[0]

    built_in_colors = [
        'b', #blue
        'g',# green
        'r',# red
        'c',# cyan
        'm',# magenta
        'y',# yellow
        'k',# black
        'w',# white
        ]


    fill_color = built_in_colors[0]

    line_styles = ['-', '--', '-.', ':', '']
    line_style = line_styles[0]

    if not ylabel:
        ylabel="y"
    if not xlabel:
        xlabel="x"
    if not label:
        label=""

    ax.plot(xs, ys, marker=symbol, zorder=1, linestyle=line_style, markerfacecolor=fill_color, label=label, markersize=marker_size, alpha=alpha, linewidth=linewidth)

    pylab.ylabel(ylabel, fontsize=font_size, fontweight=font_weight)
    pylab.xlabel(ylabel, fontsize=font_size, fontweight=font_weight)

    pylab.xticks(fontsize=font_size, fontweight=font_weight)
    pylab.yticks(fontsize=font_size, fontweight=font_weight)

    [i.set_linewidth(width) for i in ax.spines.itervalues()]
    #pylab.tick_params(which='major', length=tick_length, color=line_color, width=tick_width)
    handles, labels = ax.get_legend_handles_labels()
    ax.set_ylim((-.5, 10))
    ax.legend(handles, labels, loc = 'upper right', prop={'size':12}, numpoints = 1, frameon = False)


def lj_hard_sphere(r, sigma=1):
    T = 300
    R = 8.31e-3
    c12 = (sigma**12)*R*T
    print c12
    return c12/float(r**12)

def lj_c12_only_potential(r, c12_1, c12_2):
    return np.sqrt(c12_1*c12_2)/float(r**12)

def plot_c12_only_lj(c12=None, comp_data=None):
    c12_IPD = 3e-1 if c12 is None else c12
    c12_CD = 1e-4
    rs = np.linspace(0.3, 2, 100)
    ys_ip = [lj_c12_only_potential(r, c12_IPD, c12_IPD) for r in rs]
    print ys_ip[-1]
    ys_gr = [lj_c12_only_potential(r, c12_IPD, c12_CD) for r in rs]
    fig = create_figure()
    ax = add_axis_to_figure(fig)
    plot(ax, rs, ys_ip, label="IrPPY3-IrPPY3")
    if comp_data is not None:
        plot(ax, rs, comp_data, label="comp_data")
    plot(ax, rs, ys_gr, label="IrPPY3-Graphene")
    pylab.show()
    return ys_ip

def lj_potential(r, c12s, c6s):
    return np.sqrt(c12s[0]*c12s[1])/float(r**12) - np.sqrt(c6s[0]*c6s[1])/float(r**6)

def plot_lj_potential(ys_ip_c12_only):
    c12_IPD = 1e-1
    c6_IPD = 1e-1
    c12_CD = 1e-4
    c6_CD = 1e-4
    rs = np.linspace(0.3, 2, 100)
    ys_ip = [lj_potential(r, [c12_IPD]*2, [c6_IPD]*2) for r in rs]
    print ys_ip[-1]
    ys_gr = [lj_potential(r, [c12_IPD, c12_CD], [c6_IPD, c6_CD]) for r in rs]
    fig = create_figure()
    ax = add_axis_to_figure(fig)
    plot(ax, rs, ys_ip, label="IrPPY3-IrPPY3")
    plot(ax, rs, ys_ip_c12_only, label="IrPPY3-IrPPY3 c12 only")
    plot(ax, rs, ys_gr, label="IrPPY3-Graphene")
    pylab.show()

if __name__ == "__main__":
    ys_ip = plot_c12_only_lj(3e-1)
    #plot_c12_only_lj(4e-1, comp_data=ys_ip)
    #plot_lj_potential(ys_ip)