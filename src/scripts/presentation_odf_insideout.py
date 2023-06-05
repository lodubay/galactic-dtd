"""
Plot #2 for my AAS 241. This plot compares the [O/Fe] distributions from
APOGEE to the inside-out + power-law model at two different bins of z-height.
"""

import argparse
import matplotlib.pyplot as plt
from distribution_functions import setup_axes, plot_distributions
from ofe_distribution import vice_mdf, apogee_mdf
from utils import multioutput_to_pandas
from apogee_tools import import_apogee
import paths

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

# Plot settings
OFE_LIM = (-0.1, 0.5)
ABSZ_BINS = [0, 0.5, 2]

def main(verbose=False, cmap='plasma_r'):
    output = 'diffusion/insideout/powerlaw_slope11'
    
    # Set up plot
    fig, axs = setup_axes(ncols=2, figure_width=6, xlabel='[O/Fe]', 
                          xlim=OFE_LIM, major_tick_spacing=0.2, 
                          cmap_name=cmap, absz_bins=ABSZ_BINS,
                          panel_aspect_ratio=1.2)
    
    if verbose:
        print('Plotting VICE [O/Fe] distribution from %s...' % output)
    stars = multioutput_to_pandas(output)
    plot_distributions(vice_mdf, stars, axs[:,0], linewidth=1.5,
                       label='VICE', cmap_name=cmap, absz_bins=ABSZ_BINS)
    
    if verbose:
        print('Plotting APOGEE [O/Fe] distribution...')
    apogee_data = import_apogee()
    plot_distributions(apogee_mdf, apogee_data, axs[:,1], linewidth=1.5,
                       label='APOGEE DR17', cmap_name=cmap, absz_bins=ABSZ_BINS)
            
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'presentation_plot02.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_plot02.py',
        description='Generate [O/Fe] distributions from select VICE' + \
            ' outputs for my talk at AAS 241.'
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
        default='plasma_r',
        help='Name of colormap for color-coding by radius (default: plasma_r)')
    args = parser.parse_args()
    # main(verbose=args.verbose, overwrite=args.overwrite, cmap=args.cmap)
    main(**vars(args))
