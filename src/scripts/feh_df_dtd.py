"""
This script plots [Fe/H] distribution functions from VICE runs with the same
SFH but different DTDs.
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import distribution_functions as dfs
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee
from utils import get_color_list
import paths
import _globals

# Multizone outputs
SFH = 'insideout'
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
DTD_LABELS = ['Prompt', 
              'Powerlaw\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1$ Gyr)',
              'Triple-system']
# Plot settings
NBINS = 100
FEH_LIM = (-1.1, 0.6)
SMOOTH_WIDTH = 0.2
CMAP = 'plasma_r'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    apogee_data = import_apogee()
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(DTD_LIST)+1, 
                              figure_width=_globals.TWO_COLUMN_WIDTH, 
                              cmap_name=CMAP, xlabel='[Fe/H]', xlim=FEH_LIM, 
                              major_tick_spacing=0.5)
    fig.subplots_adjust(top=0.9, left=0.07, right=0.98, bottom=0.25)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    # plot
    mdf_kwargs = {'bins': NBINS, 'range': FEH_LIM, 'smoothing': SMOOTH_WIDTH}
    for i, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH, dtd, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[fe/h]', colors, 
                                label=DTD_LABELS[i], **mdf_kwargs)
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'FE_H', colors, **mdf_kwargs)
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'feh_df_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
