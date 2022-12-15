"""
Plot #1 for my AAS 241 talk. This script plots [O/Fe] vs [Fe/H] of stellar
particles from a VICE multi-zone run with the inside-out SFH and power law DTD.
Two regions are included: midplane and out-of-plane, both in the solar annulus.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from utils import import_allStar, multioutput_to_pandas, \
    filter_multioutput_stars, sample_dataframe
from ofe_feh_apogee import plot_contours
import paths
from _globals import ZONE_WIDTH

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

# Plot parameters
FEH_LIM = (-1.3, 0.6)
OFE_LIM = (-0.25, 0.65)
ABSZ_BINS = [(0.5, 2), (0, 0.5)] # kpc
GALR_LIM = (7, 9) # kpc

def main(verbose=False, overwrite=False, cmap='winter'):
    # APOGEE data
    if verbose:
        print('Importing APOGEE allStar data...')
    apogee_data = import_allStar(only_giants=True)
    
    # VICE output
    if verbose:
        print('Importing VICE stellar population data...')
    stars = multioutput_to_pandas('diffusion/insideout/powerlaw_slope11')
    
    # Set up figure and axes
    fig, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True, sharey=True)
    plt.subplots_adjust(left=0.13, bottom=0.09, top=0.95, right=0.82, hspace=0.)
    
    # Set up colorbar
    subset = filter_multioutput_stars(stars, galr_lim=GALR_LIM, absz_lim=(0, 2),
                                      zone_width=ZONE_WIDTH)
    norm = Normalize(vmin=subset['zone_origin'].min() * ZONE_WIDTH, 
                     vmax=subset['zone_origin'].max() * ZONE_WIDTH)
    cax = plt.axes([0.85, 0.09, 0.04, 0.86])
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label(r'Galactocentric radius at birth [kpc]', labelpad=5)
    cax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    galr_lim = GALR_LIM
    
    for ax, absz_lim in zip(axs, ABSZ_BINS):
        # Scatter plot of VICE abundances
        subset = filter_multioutput_stars(stars, galr_lim, absz_lim, ZONE_WIDTH)
        # weight random sample by particle mass
        sample_weights = subset['mass'] / subset['mass'].sum()
        sample = sample_dataframe(subset, 10000, weights=sample_weights)
        # Scatter plot of random sample of stellar particles
        scatters = ax.scatter(sample['[fe/h]'], sample['[o/fe]'], s=0.5,
                              c=sample['zone_origin'] * ZONE_WIDTH, cmap=cmap,
                              # norm=norm, 
                              rasterized=True, edgecolor='none')
        
        # APOGEE contours
        contours = plot_contours(ax, apogee_data, overwrite=overwrite,
                                 absz_lim=absz_lim, galr_lim=galr_lim, 
                                 colors='k', linestyles=[':', '--', '-'],
                                 linewidths=0.5)
        
        # Label z-height bins
        ax.text(0.95, 0.92, r'%s kpc $\leq |z| <$ %s kpc' % absz_lim,
                transform=ax.transAxes, va='top', ha='right')
        
    # Contour legend
    contour_labels = [r'$1\sigma$', r'$2\sigma$', r'$3\sigma$']
    # for i, label in enumerate(contour_labels):
    #     contours.collections[i].set_label(label)
    axs[1].legend(contours.legend_elements()[0], contour_labels, 
                  loc='lower left', title='APOGEE\ncontours', frameon=False)
        
    # Configure axes
    axs[0].set_title(r'%s kpc $\leq R_{\rm{Gal}} <$ %s kpc' % GALR_LIM)
    
    axs[1].set_xlabel('[Fe/H]', labelpad=5)
    axs[0].set_xlim(FEH_LIM)
    axs[0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0].xaxis.set_minor_locator(MultipleLocator(0.1))
    
    axs[0].set_ylabel('[O/Fe]')
    axs[1].set_ylabel('[O/Fe]')
    axs[0].set_ylim(OFE_LIM)
    axs[0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    plt.savefig(paths.figures / 'presentation_plot01.png', dpi=300)


def normalize_colorbar(stars):
    """
    Calculate the colorbar normalization for the full migration output data
    before subsampling.

    Parameters
    ----------
    stars : pandas DataFrame
        VICE migration output

    Returns
    -------
    norm : instance of matplotlib.colors.Normalize
    """
    # Filter to full range to be plotted
    return norm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_plot01.py',
        description='Generate a plot of [O/Fe] vs [Fe/H] from select VICE' + \
            ' outputs for my talk at AAS 241.'
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-o', '--overwrite', action='store_true',
        help='Whether to overwrite saved 2D KDE data (takes longer)')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
        default='winter',
        help='Name of colormap for color-coding VICE output (default: winter)')
    args = parser.parse_args()
    # main(verbose=args.verbose, overwrite=args.overwrite, cmap=args.cmap)
    main(**vars(args))
