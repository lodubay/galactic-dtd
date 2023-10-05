"""
This script plots [O/Fe] distribution functions from VICE runs with the same
SFH but different DTDs.
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import distribution_functions as dfs
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee
from utils import get_color_list, highlight_panels
import paths
import _globals

# Multizone outputs
SFH = 'earlyburst'
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
DTD_LABELS = ['Two-population', 
              'Power-law\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1$ Gyr)',
              'Triple-system']
# Plot settings
NBINS = 100
OFE_LIM = (-0.15, 0.55)
SMOOTH_WIDTH = 0.05
CMAP = 'plasma_r'

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    apogee_data = import_apogee()
    # Limit size of plot in poster format
    if style == 'poster':
        dtd_list = ['powerlaw_slope11', 'exponential_timescale15', 'plateau_width10']
        figwidth = _globals.ONE_COLUMN_WIDTH * 1.8
    else:
        dtd_list = DTD_LIST
        figwidth = _globals.TWO_COLUMN_WIDTH
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(dtd_list)+1, 
                              figure_width=figwidth, 
                              cmap_name=CMAP, xlabel='[O/Fe]', xlim=OFE_LIM, 
                              major_tick_spacing=0.2, major_minor_tick_ratio=4.,
                              cbar_width=0.4)
    fig.subplots_adjust(top=0.9, left=0.07, right=0.98, bottom=0.25)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    # plot
    mdf_kwargs = {'bins': NBINS, 'range': OFE_LIM, 'smoothing': SMOOTH_WIDTH}
    for i, dtd in enumerate(dtd_list):
        output_name = '/'.join(['gaussian', SFH, dtd, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[o/fe]', colors, 
                                label=DTD_LABELS[i], **mdf_kwargs)
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'O_FE', colors, **mdf_kwargs)
    highlight_panels(fig, axs, [(0,-1),(1,-1),(2,-1)])
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'ofe_df_dtd')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_df_dtd.py',
        description='Compare distribution functions of [O/Fe] from many ' +
        'Galactic regions for VICE outputs with different delay time ' +
        'distributions.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
