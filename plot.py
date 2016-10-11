from os.path import dirname
import os
import matplotlib
if not os.environ.has_key("DISPLAY"):
    matplotlib.use("Agg")
import pylab
#matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='Arial')
matplotlib.rc('text', usetex='false')

LINE_STYLES = ['-', '--', '-.', ':', '']

def in_dict(key, dictionary):
    if dictionary:
        return key in dictionary and dictionary[key]
    return False

def get_nice_colours(color_name):
    color_dict = {
        'g': '#32CD32',
        'y': '#ffd700',
        'o': '#ffa500',
    }
    if not in_dict(color_name, color_dict):
        return color_name
    return color_dict[color_name]


def create_figure(figsize=(10,8)):
    fig = pylab.figure(figsize=figsize)
    fig.hold = True
    return fig

def create_subplot_layout(number_of_plots, current_plot, dense=True):
    rows = 1
    cols = 1
    while rows*cols < number_of_plots:
        rows += 1
        if rows*cols < number_of_plots:
            cols += 1
    if dense:
        return int('{0}{1}{2}'.format(rows, cols, current_plot))
    return rows, cols, current_plot


def add_axis_to_figure(fig, subplot_layout=111, sharex=None, sharey=None):
    return fig.add_subplot(subplot_layout, sharex=sharex, sharey=sharey)

def save_figure(fig, fig_name, image_format='png', dpi=300):
    graph_dir = dirname(fig_name)
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)
    file_name = '{0}.{1}'.format(fig_name, image_format)
    #print "saving to file: {0}".format(file_name)
    fig.savefig(file_name, dpi=dpi)
    pylab.close()


def plot_batch(ax, data):
    for label, (xs, ys, color) in data.items():
        plot(ax, xs, ys, line_color=color, fill_color=color, label=label)

def plot(ax,
         xs,
         ys,
         color='black',
         title="",
         label=None,
         xlabel="",
         ylabel="",
         zorder=None,
         alpha=1.0,
         symbol=None,
         linewidth=0.5,
         xlim=None,
         ylim=None,
         line_color=None,
         fill_color=None,
         line_style=LINE_STYLES[0],
         marker_size=5,
         marker_outline_color=None,
         drawstyle=None, #'default',
         legend_position=None, #'upper right',
         legend_frame=False,
         font_size=12,
         axis_line_width=1,
         tick_width=0.75,
         tick_length=4,
         font_weight="normal",
         dashes=None,
         text=None,
         text_location=(0.82, 0.9),
         ):

    if line_color is None:
        line_color = color
    if fill_color is None:
        fill_color = color
    if marker_outline_color is None:
        marker_outline_color = color
    additional_args = {}
    if dashes is not None:
        additional_args["dashes"] = dashes
    ax.plot(xs, ys, marker=symbol, zorder=zorder, linestyle=line_style, color=get_nice_colours(line_color),
            markerfacecolor=get_nice_colours(fill_color), label=label, markersize=marker_size, linewidth=linewidth,
            alpha=alpha, markeredgecolor=get_nice_colours(marker_outline_color), drawstyle=drawstyle, **additional_args)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if title:
        ax.set_title(title, fontsize=font_size+2)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=font_size, fontweight=font_weight)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=font_size, fontweight=font_weight)

    ax.tick_params(axis='both', which='both', labelsize=font_size-2)
    #ax.tick_params(axis='both', which='minor', labelsize=8)

    [i.set_linewidth(axis_line_width) for i in ax.spines.itervalues()]
    ax.tick_params(which='major', length=tick_length, color="k", width=tick_width)
    handlers, labels = remove_label_duplicates(*ax.get_legend_handles_labels())
    ax.legend(handlers, labels, loc=legend_position, prop={'size': font_size-1}, numpoints=1, frameon=legend_frame)
    if text:
        ax.text(text_location[0], text_location[1], text, transform=ax.transAxes, fontsize=font_size, horizontalalignment="center", verticalalignment="center")


def remove_label_duplicates(handlers, labels):
    updated_data = {}
    if handlers and labels:
        for l, h in zip(labels, handlers):
            if l in updated_data:
                if h._alpha > updated_data[l]._alpha:
                    updated_data[l] = h
            else:
                updated_data[l] = h
        labels, handlers = zip(*updated_data.items())
    return handlers, labels

if __name__=="__main__":
    fig = create_figure((8,7))
    ax = add_axis_to_figure(fig)
    plot(ax, [0,1], [0,1], xlabel="x", ylabel="y", label="test", line_color="g", symbol="o", linewidth=2, axis_line_width=1.5, tick_width=1.2, font_size=14)
    pylab.show()