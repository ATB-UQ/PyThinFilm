import matplotlib.pyplot as p
from PIL import Image
import numpy as np
from os.path import basename
import glob
import matplotlib
fontsize = 11
name_map = {"n15": "2.0", "n40": "5.3", "n108": "14.1", "n209":"26.3"}
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='Times New Roman')
matplotlib.rc('text', usetex='false')
def figure_collage(png_files):

    subplot_dim = int(np.ceil(np.sqrt(len(png_files))))

    fig, axarr = p.subplots(*[subplot_dim]*2, figsize=(8.2, 8.8 ))

    def indices_for_fig(n):
        return ((n // subplot_dim), n - (n // subplot_dim) * subplot_dim)
    show_y_axis_label = [0,4,8,12]
    show_title = [0,1,2,3]
    for (n, png_file) in enumerate(sorted(png_files, key=lambda x:float(basename(x)[:-4].split("_")[1][1:]))):
        image = Image.open(png_file)
        image = image.crop((405, 24, 1850, 1680))
        axarr[indices_for_fig(n)].imshow(image)
        image_words = basename(png_file[:-4]).split("_")
        wtp, cutoff = name_map[image_words[1]], image_words[2]
        if n in show_title:
            axarr[indices_for_fig(n)].set_title(
                '{0} (nm) cutoff'.format(cutoff),
                fontsize=fontsize,
            )
        axarr[indices_for_fig(n)].set_axis_off()
        if n in show_y_axis_label:
            axarr[indices_for_fig(n)].text(-.2, 0.5, "{0} wt%".format(wtp), transform=axarr[indices_for_fig(n)].transAxes, fontsize=fontsize, horizontalalignment="center", verticalalignment="center")
    for n in range(len(png_files), subplot_dim**2):
        axarr[indices_for_fig(n)].set_axis_off()

    p.tight_layout()
    fig.subplots_adjust(wspace=.05, hspace=.05, right=.99, left=0.09, top=0.97, bottom=0.01)

    #p.show()
    fig.savefig('collage.png', dpi=600)

figure_collage(glob.glob("images/whole_*.png"))