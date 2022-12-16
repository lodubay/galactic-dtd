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
FEH_LIM = (-1.4, 0.6)
OFE_LIM = (-0.25, 0.65)
ABSZ_BINS = [(0, 0.5), (0.5, 2)] # kpc
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
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    bounds = {
        'left': 0.1,
        'bottom': 0.18,
        'right': 0.85,
        'top': 0.88,
        'hspace': 0.05,
        'wspace': 0.05
    }
    plt.subplots_adjust(**bounds)
    
    # Set up colorbar
    # subset = filter_multioutput_stars(stars, galr_lim=GALR_LIM, absz_lim=(0, 2),
    #                                   zone_width=ZONE_WIDTH)
    norm = Normalize(vmin=0, vmax=16)
    cax = plt.axes([bounds['right'] + bounds['wspace']/2.5, bounds['bottom'], 
                    0.03, bounds['top'] - bounds['bottom']])
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label(r'Birth radius [kpc]', labelpad=6)
    cax.yaxis.set_major_locator(MultipleLocator(4))
    cax.yaxis.set_minor_locator(MultipleLocator(1))
    
    galr_lim = GALR_LIM
    
    for ax, absz_lim in zip(axs, ABSZ_BINS):
        # Scatter plot of VICE abundances
        subset = filter_multioutput_stars(stars, galr_lim, absz_lim, ZONE_WIDTH)
        # weight random sample by particle mass
        sample_weights = subset['mass'] / subset['mass'].sum()
        sample = sample_dataframe(subset, 10000, weights=sample_weights)
        # Scatter plot of random sample of stellar particles
        scatters = ax.scatter(sample['[fe/h]'], sample['[o/fe]'], s=1,
                              c=sample['zone_origin'] * ZONE_WIDTH, cmap=cmap,
                              # norm=norm, 
                              rasterized=True, edgecolor='none')
        
        # APOGEE contours
        contours = plot_contours(ax, apogee_data, overwrite=overwrite,
                                 absz_lim=absz_lim, galr_lim=galr_lim, 
                                 colors='k', linestyles=[':', '--', '-'],
                                 linewidths=0.5)
        
        # Label z-height bins
        ax.set_title(r'%s kpc $\leq |z| <$ %s kpc' % absz_lim, pad=12)
        
    # Contour legend
    contour_labels = [r'$1\sigma$', r'$2\sigma$', r'$3\sigma$']
    # for i, label in enumerate(contour_labels):
    #     contours.collections[i].set_label(label)
    contours.legend_elements()[0].reverse()
    axs[0].legend(contours.legend_elements()[0], contour_labels, 
                  loc='lower left', title='APOGEE\ncontours', frameon=False)
        
    # Configure axes
    axs[1].text(0.95, 0.92, r'%s kpc $\leq R_{\rm{Gal}} <$ %s kpc' % GALR_LIM,
                transform=ax.transAxes, va='top', ha='right')
    # axs[0].set_title(r'%s kpc $\leq R_{\rm{Gal}} <$ %s kpc' % GALR_LIM)
    
    axs[0].set_xlabel('[Fe/H]', labelpad=6)
    axs[1].set_xlabel('[Fe/H]', labelpad=6)
    axs[0].set_xlim(FEH_LIM)
    axs[0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0].xaxis.set_minor_locator(MultipleLocator(0.1))
    
    axs[0].set_ylabel('[O/Fe]', labelpad=-6)
    axs[0].set_ylim(OFE_LIM)
    axs[0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    plt.savefig(paths.figures / 'presentation_ofe_feh_insideout.png', dpi=300)
    plt.close()


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
