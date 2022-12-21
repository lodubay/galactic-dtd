"""
Plot #4 for my AAS 241. This plot compares the [Fe/H] distributions from
APOGEE to the Conroy22 + power-law model at two different bins of z-height.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
from distribution_functions import setup_axes, plot_distributions
from feh_distribution import vice_mdf, apogee_mdf
from utils import import_allStar, multioutput_to_pandas, discrete_colormap
from _globals import GALR_BINS
import paths

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

# Plot settings
FEH_LIM = (-1.3, 0.8)
ABSZ_BINS = [0, 0.5, 2]

def main(verbose=False, cmap_name='viridis_r'):
    output = 'diffusion/conroy22/powerlaw_slope11'
    
    # Set up plot
    fig, axs = plt.subplots(len(ABSZ_BINS)-1, 2, sharex=True, sharey='row',
                            figsize=(9, 4.5))
    bounds = {
        'left': 0.08,
        'bottom': 0.15,
        'right': 0.88,
        'top': 0.92,
        'hspace': 0.05,
        'wspace': 0.05
    }
    plt.subplots_adjust(**bounds)
    
    # Format x-axis
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    
    # Remove spines and y-axis labels
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
        # Set bottom ticks pointing out
        ax.tick_params(axis='x', which='both', direction='out')
    # Add common y-axis label
    fig.text(0.01, 0.5, r'Distance from midplane $|z|$',
             ha='left', va='center', rotation='vertical', 
             size=plt.rcParams['axes.labelsize'])
    # Label rows
    for i in range(len(ABSZ_BINS)-1):
        absz_lim = tuple(ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i])
        axs[i,0].set_ylabel(r'$%s - %s$ kpc' % absz_lim, 
                            size=plt.rcParams['ytick.labelsize'])
    
    # Add colorbar
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = plt.axes([bounds['right'] + bounds['wspace']/2.5, bounds['bottom'], 
                    0.02, bounds['top'] - bounds['bottom']])
    # cax = plt.axes([0.5 - (cbar_width / 2), 0.09, cbar_width, 0.02])
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label('Galactocentric radius [kpc]')

    if verbose:
        print('Plotting APOGEE [Fe/H] distribution...')
    apogee_data = import_allStar()
    plot_distributions(apogee_mdf, apogee_data, axs[:,0], linewidth=1.5,
                       label='APOGEE', cmap_name=cmap, absz_bins=ABSZ_BINS,
                       func_kwargs={'xlim': FEH_LIM})
    
    if verbose:
        print('Plotting VICE [Fe/H] distribution from %s...' % output)
    stars = multioutput_to_pandas(output)
    plot_distributions(vice_mdf, stars, axs[:,1], linewidth=1.5,
                       label=r'Power law ($t^{-1.1}$)', cmap_name=cmap, 
                       absz_bins=ABSZ_BINS, func_kwargs={'xlim': FEH_LIM})
            
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / 'presentation_mdf_conroy22.png', dpi=300)
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
    main(verbose=args.verbose, cmap_name=args.cmap)
