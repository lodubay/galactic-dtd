"""
This script compares the MDFs between the analog and gaussian migration schemes.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
from multizone_stars import MultizoneStars
import distribution_functions as dfs
from apogee_tools import import_apogee
from utils import box_smooth, get_color_list, discrete_colormap
import paths
import _globals

MIGRATION_SCHEMES = ['diffusion', 'gaussian']
ROW_LABELS = ['h277 analogue', 'Gaussian sampling']
# Plot settings
NBINS = 100
FEH_LIM = (-1.2, 0.6)
SMOOTH_WIDTH = 0.2
CMAP = 'plasma_r'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    apogee_data = import_apogee()
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(MIGRATION_SCHEMES)+1, 
                              figure_width=_globals.TWO_COLUMN_WIDTH, 
                              cmap=CMAP, xlabel='[Fe/H]', xlim=FEH_LIM, 
                              major_tick_spacing=0.5, cbar_width=0.6,
                              include_yaxis=True)
    # fig.subplots_adjust(top=0.92, left=0.04, right=0.96)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    mdf_kwargs = {'bins': NBINS, 'range': FEH_LIM, 'smoothing': SMOOTH_WIDTH}
    # plot varying SFH, constant DTD
    # title_size = plt.rcParams['axes.titlesize']
    # fig.text(0.22, 0.92, r'(Exponential DTD, $\tau_{\rm Ia}=1.5$ Gyr)',
    #          ha='center', va='bottom')
    for i, scheme in enumerate(MIGRATION_SCHEMES):
        output_name = '/'.join((scheme, 'insideout/powerlaw_slope11/diskmodel'))
        mzs = MultizoneStars.from_output(output_name)
        # mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[fe/h]', colors, #titlepad=-18,
                                label=ROW_LABELS[i], **mdf_kwargs)
    # right column APOGEE for comparison
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'FE_H', colors, #titlepad=0,
                          **mdf_kwargs)
    
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'migration_mdf.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()
