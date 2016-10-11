from plot import create_figure, add_axis_to_figure, plot, save_figure
import numpy as np
DEPOSITED_DATA_PATH = "rdf_xvg_data/n{0}/rdf_IP_IP.xvg"
RANDOM_DATA_PATH = "rdf_xvg_data/n{0}/rdf_IPD_IPD.xvg"
TITLE_MAP = {15:"2.0 wt%",40:"5.3 wt%",108:"14.1 wt%",209:"26.3 wt%"}

def smooth(x,window_len=20,window_type='flat'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window_type in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window_type == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window_type+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def plot_single_rdf(deposited, random, n_irppy):
    fig = create_figure((4,3))
    ax = add_axis_to_figure(fig)
    plot_to_axis(ax, deposited, random, title=TITLE_MAP[n_irppy])
    fig.tight_layout()
    save_figure(fig, "./rdf_{0}".format(TITLE_MAP[n_irppy].replace(" ", "_")), image_format="eps")

def plot_to_axis(ax, deposited, random, hide_yaxis_label=False, hide_xaxis_label=False, title=None):
    ylabel = None if hide_yaxis_label else "g(r)"
    xlabel = None if hide_xaxis_label else "r (nm)"
    #plot(ax, deposited[0], deposited[1], color="k", label="Ir(ppy)3", xlabel="r (nm)", ylabel=ylabel, zorder=2, linewidth=1, title=title, legend_frame=False)
    #plot(ax,  random[0], random[1], color="k", label="Random", zorder=1, linewidth=1, xlim=(0,4), legend_position="upper right", dashes=(4,2), legend_frame=False)
    plot(ax, deposited[0], deposited[1], color="b", xlabel=xlabel, ylabel=ylabel, zorder=2, linewidth=1.5, text=title, legend_frame=False)
    plot(ax,  random[0], random[1], color="k", zorder=1, linewidth=1, xlim=(0,4), legend_position="upper right", dashes=(4,2), legend_frame=False)

def plot_all_on_single_fig():
    fig = create_figure((6,5))
    ax = add_axis_to_figure(fig, subplot_layout=int("221"))
    deposited = parse_xvg(DEPOSITED_DATA_PATH.format(15))
    random = parse_xvg(RANDOM_DATA_PATH.format(15), smooth_data=True)
    plot_to_axis(ax, deposited, random, title=TITLE_MAP[15], hide_xaxis_label=True)
    for i, n_irppy in enumerate((40, 108, 209)):
        deposited = parse_xvg(DEPOSITED_DATA_PATH.format(n_irppy))
        random = parse_xvg(RANDOM_DATA_PATH.format(n_irppy))
        ax = add_axis_to_figure(fig, subplot_layout=int("22{0}".format(i+2)), sharey=ax)
        hide_y_label = i in [0, 2]
        hide_x_label = i in [0]
        plot_to_axis(ax, deposited, random, hide_yaxis_label=hide_y_label, hide_xaxis_label=hide_x_label, title=TITLE_MAP[n_irppy])
    fig.tight_layout()
    save_figure(fig, "./rdf_combined", image_format="eps")

def parse_xvg(data_path, smooth_data=False):
    with open(data_path) as fh:
        xvg_str = fh.read()
    xs, ys = zip(*[map(float, l.split()) for l in xvg_str.splitlines() if l[0] not in ["#", "@"]])
    if smooth_data:
        tail_cutoff = 1.1
        y_tail = [y for x,y in zip(xs, ys) if x >= tail_cutoff]
        y_head = [y for x,y in zip(xs, ys) if x < tail_cutoff]
        y_tail = smooth(np.array(y_tail), 8, "flat")
        ys = y_head + list(y_tail)
    return xs, ys

def plot_all():
    for n_irppy in (15, 40, 108, 209):
        deposited = parse_xvg(DEPOSITED_DATA_PATH.format(n_irppy))
        random = parse_xvg(RANDOM_DATA_PATH.format(n_irppy), smooth_data=True)
        plot_single_rdf(deposited, random, n_irppy)

if __name__=="__main__":
    #plot_all_on_single_fig()
    plot_all()
