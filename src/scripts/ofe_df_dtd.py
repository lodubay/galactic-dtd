"""
This script plots [O/Fe] distribution functions from VICE runs with the same
SFH_DEFAULT but different DTDs.
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
SFH_DEFAULT = 'earlyburst'
SFH_LABEL = 'Early-burst SFH'
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

def main(style='paper', sfh=SFH_DEFAULT):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    apogee_data = import_apogee()
    # Limit size of plot in poster format
    if style == 'poster':
        dtd_list = DTD_LIST[1:4]
        dtd_labels = DTD_LABELS[1:4]
        figwidth = _globals.ONE_COLUMN_WIDTH * 1.8
        cbar_width = 0.4
    else:
        dtd_list = DTD_LIST
        dtd_labels = DTD_LABELS
        figwidth = _globals.TWO_COLUMN_WIDTH
        cbar_width = 0.4
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(dtd_list)+1, 
                              figure_width=figwidth, 
                              cmap=CMAP, xlabel='[O/Fe]', xlim=OFE_LIM, 
                              major_tick_spacing=0.2, major_minor_tick_ratio=4.,
                              cbar_width=cbar_width, panel_aspect_ratio=1.3)
    fig.subplots_adjust(top=0.85, left=0.04, right=0.96, bottom=0.23)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    # plot
    mdf_kwargs = {'bins': NBINS, 'range': OFE_LIM, 'smoothing': SMOOTH_WIDTH}
    for i, dtd in enumerate(dtd_list):
        output_name = '/'.join(['gaussian', sfh, dtd, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[o/fe]', colors, 
                                label=dtd_labels[i], **mdf_kwargs)
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'O_FE', colors, **mdf_kwargs)
    highlight_panels(fig, axs, [(0,-1),(1,-1),(2,-1)])
    # Figure title indicating SFH model used
    fig.suptitle(SFH_LABEL)
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    fname = 'ofe_df_dtd'
    if sfh != SFH_DEFAULT:
        fname += f'_{sfh}'
    plt.savefig(paths.figures / fname)
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
    parser.add_argument('--sfh', 
                        choices=['insideout', 'lateburst', 'earlyburst', 'twoinfall'],
                        default=SFH_DEFAULT,
                        help='Which SFH model to plot for all panels.')
    args = parser.parse_args()
    main(**vars(args))
