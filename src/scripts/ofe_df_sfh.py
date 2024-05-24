"""
This script plots [O/Fe] distribution functions from VICE runs with the same
DTD but different SFHs.
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import distribution_functions as dfs
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee
from utils import get_color_list, highlight_panels
from colormaps import paultol
import paths
import _globals

# Multizone outputs
SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
SFH_LABELS = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD = 'exponential_timescale15'
DTD_LABEL = 'Exponential DTD ($\\tau=1.5$ Gyr)'
# Plot settings
NBINS = 100
OFE_LIM = (-0.15, 0.55)
SMOOTH_WIDTH = 0.05
CMAP = plt.get_cmap('plasma_r')
# CMAP = paultol.ylorbr_short

def main(style='poster'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    apogee_data = import_apogee()
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(SFH_LIST)+1, 
                              figure_width=_globals.TWO_COLUMN_WIDTH, 
                              cmap=CMAP, xlabel='[O/Fe]', xlim=OFE_LIM, 
                              major_tick_spacing=0.2, major_minor_tick_ratio=4.,
                              cbar_width=0.4, panel_aspect_ratio=1.4)
    fig.subplots_adjust(top=0.86, left=0.04, right=0.96, bottom=0.23)
    colors = get_color_list(CMAP, _globals.GALR_BINS)
    # plot
    mdf_kwargs = {'bins': NBINS, 'range': OFE_LIM, 'smoothing': SMOOTH_WIDTH}
    for i, sfh in enumerate(SFH_LIST):
        output_name = '/'.join(['gaussian', sfh, DTD, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[o/fe]', colors, 
                                label=SFH_LABELS[i], **mdf_kwargs)
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'O_FE', colors, **mdf_kwargs)
    highlight_panels(fig, axs, [(0,-1),(1,-1),(2,-1)])
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    fig.suptitle(DTD_LABEL)
    plt.savefig(paths.figures / 'ofe_df_sfh', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_df_sfh.py',
        description='Compare distribution functions of [O/Fe] from many ' +
        'Galactic regions for VICE outputs with different SFHs.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
