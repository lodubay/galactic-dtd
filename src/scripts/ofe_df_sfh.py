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
from utils import get_color_list
import paths
import _globals

# Multizone outputs
SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
SFH_LABELS = ['Inside-Out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD = 'exponential_timescale15'
# Plot settings
NBINS = 100
OFE_LIM = (-0.15, 0.55)
SMOOTH_WIDTH = 0.05
CMAP = 'plasma_r'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    apogee_data = import_apogee()
    # Set up plot
    fig, axs = dfs.setup_axes(ncols=len(SFH_LIST)+1, 
                              figure_width=_globals.TWO_COLUMN_WIDTH, 
                              cmap_name=CMAP, xlabel='[O/Fe]', xlim=OFE_LIM, 
                              major_tick_spacing=0.2, major_minor_tick_ratio=4.,
                              cbar_width=0.4)
    fig.subplots_adjust(top=0.9, left=0.07, right=0.98)
    colors = get_color_list(plt.get_cmap(CMAP), _globals.GALR_BINS)
    # plot
    mdf_kwargs = {'bins': NBINS, 'range': OFE_LIM, 'smoothing': SMOOTH_WIDTH}
    for i, sfh in enumerate(SFH_LIST):
        output_name = '/'.join(['gaussian', sfh, DTD, 'diskmodel'])
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data, inplace=True)
        dfs.plot_multizone_mdfs(mzs, axs[:,i], '[o/fe]', colors, 
                                label=SFH_LABELS[i], **mdf_kwargs)
    dfs.plot_apogee_mdfs(apogee_data, axs[:,-1], 'O_FE', colors, **mdf_kwargs)
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'ofe_df_sfh.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
