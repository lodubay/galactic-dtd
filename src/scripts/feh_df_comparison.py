"""
This script plots [Fe/H] distribution functions from VICE runs with the same
DTD but different SFHs.
"""

import numpy as np
import matplotlib.pyplot as plt
import distribution_functions as dfs
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee
from utils import get_color_list, highlight_panels
import paths
import _globals

# Multizone outputs
SFH_LIST = ['insideout', 'twoinfall']
SFH_LABELS = ['Inside-out SFH', 'Two-infall SFH']
DTD_MODEL = 'exponential_timescale15' # hold constant while varying SFH
DTD_LIST = ['powerlaw_slope14',
            'exponential_timescale30']
DTD_LABELS = ['Power-law DTD\n($\\alpha=-1.4$)',
              'Exponential DTD\n($\\tau=3$ Gyr)']
# DTD_LIST = ['prompt', 'triple']
# DTD_LABELS = ['Two-population DTD', 'Triple-system DTD']
SFH_MODEL = 'insideout' # hold constant while varying DTD
# Plot settings
NBINS = 100
FEH_LIM = (-1.2, 0.6)
SMOOTH_WIDTH = 0.2
CMAP = 'plasma_r'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    apogee_data = import_apogee()
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(SFH_LIST)+len(DTD_LIST)+1, 
                              figure_width=_globals.TWO_COLUMN_WIDTH, 
                              cmap=CMAP, xlabel='[Fe/H]', xlim=FEH_LIM, 
                              major_tick_spacing=0.5, cbar_width=0.4)
    fig.subplots_adjust(top=0.92, left=0.04, right=0.96)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    mdf_kwargs = {'bins': NBINS, 'range': FEH_LIM, 'smoothing': SMOOTH_WIDTH}
    # plot varying SFH, constant DTD
    # title_size = plt.rcParams['axes.titlesize']
    # fig.text(0.22, 0.92, r'(Exponential DTD, $\tau_{\rm Ia}=1.5$ Gyr)',
    #          ha='center', va='bottom')
    for i, sfh in enumerate(SFH_LIST):
        output_name = '/'.join(['gaussian', sfh, DTD_MODEL, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[fe/h]', colors, #titlepad=-18,
                                label=SFH_LABELS[i], **mdf_kwargs)
    # center column APOGEE for comparison
    dfs.plot_apogee_mdfs(apogee_data, axs[:,2], 'FE_H', colors, #titlepad=0,
                          **mdf_kwargs)
    # plot varying DTD, constant SFH
    # fig.text(0.78, 0.92, '(Inside-out SFH)', ha='center', va='bottom')
    for i, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH_MODEL, dtd, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,3+i], '[fe/h]', colors, #titlepad=-12,
                                label=DTD_LABELS[i], **mdf_kwargs)
    highlight_panels(fig, axs, [(0,2),(1,2),(2,2)])
    
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'feh_df_comparison.pdf', dpi=300)
    plt.close()


def get_subplot_edges(fig, axs):
    # Get the bounding boxes of the axes including text decorations
    r = fig.canvas.get_renderer()
    get_bbox = lambda ax: ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
    bboxes = list(map(get_bbox, axs.flat))
    
    #Get the minimum and maximum extent, get the coordinate half-way between those
    xmax = np.array(list(map(lambda b: b.x1, bboxes))).reshape(axs.shape).max(axis=0)
    xmin = np.array(list(map(lambda b: b.x0, bboxes))).reshape(axs.shape).min(axis=0)
    xs = np.c_[xmax[1:], xmin[:-1]].mean(axis=1)
    return xs


if __name__ == '__main__':
    main()
