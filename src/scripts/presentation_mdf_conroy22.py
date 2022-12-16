"""
Plot #4 for my AAS 241. This plot compares the [Fe/H] distributions from
APOGEE to the Conroy22 + power-law model at two different bins of z-height.
"""

import argparse
import matplotlib.pyplot as plt
from distribution_functions import setup_axes, plot_distributions
from feh_distribution import vice_mdf, apogee_mdf
from utils import import_allStar, multioutput_to_pandas
import paths

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

# Plot settings
FEH_LIM = (-1.1, 0.6)
ABSZ_BINS = [0, 0.5, 2]

def main(verbose=False, cmap='viridis_r'):
    output = 'diffusion/conroy22/powerlaw_slope11'
    
    # Set up plot
    fig, axs = setup_axes(ncols=2, figure_width=6, xlabel='[Fe/H]', 
                          xlim=FEH_LIM, major_tick_spacing=0.5, 
                          cmap_name=cmap, absz_bins=ABSZ_BINS,
                          panel_aspect_ratio=1.3)
    fig.subplots_adjust(left=0.12)
    
    if verbose:
        print('Plotting VICE [Fe/H] distribution from %s...' % output)
    stars = multioutput_to_pandas(output)
    plot_distributions(vice_mdf, stars, axs[:,0], linewidth=1.5,
                       label='VICE', cmap_name=cmap, absz_bins=ABSZ_BINS)
    
    if verbose:
        print('Plotting APOGEE [Fe/H] distribution...')
    apogee_data = import_allStar()
    plot_distributions(apogee_mdf, apogee_data, axs[:,1], linewidth=1.5,
                       label='APOGEE DR17', cmap_name=cmap, absz_bins=ABSZ_BINS)
            
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'presentation_plot04.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_plot04.py',
        description='Generate [Fe/H] distributions from select VICE' + \
            ' outputs for my talk at AAS 241.'
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
        default='viridis_r',
        help='Name of colormap for color-coding by radius (default: viridis_r)')
    args = parser.parse_args()
    # main(verbose=args.verbose, overwrite=args.overwrite, cmap=args.cmap)
    main(**vars(args))
