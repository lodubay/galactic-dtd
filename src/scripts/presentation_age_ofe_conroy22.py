"""
Plot for my AAS 241 talk. This script plots [O/Fe] vs the age of stellar
particles from a VICE multi-zone run with the Conroy+ 2022 SFH and two
different DTDs.
Two regions are included: midplane and out-of-plane, both in the solar annulus.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from utils import import_astroNN, apogee_region, multioutput_to_pandas, \
    filter_multioutput_stars
from scatter_plot_grid import setup_colorbar, plot_vice_sample
from age_ofe import plot_vice_medians, plot_astroNN_medians
import paths
from _globals import ZONE_WIDTH

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

OUTPUTS = ['diffusion/conroy22/powerlaw_slope11',
           # 'diffusion/conroy22/exponential_timescale30',
            'diffusion/conroy22/plateau_width300_slope11',
           ]
LABELS = [r'Power law ($t^{-1.1}$)',
          # r'Exponential ($\tau=3$ Gyr)',
          r'Plateau (300 Myr)',
          ]

# Plot parameters
AGE_LIM_LINEAR = (-2, 12)
AGE_LIM_LOG = (0.2, 20)
OFE_LIM = (-0.15, 0.65)
ABSZ_BINS = [(0.5, 2), (0, 0.5)] # kpc
GALR_LIM = (7, 9) # kpc

def main(verbose=False, log=False, cmap='winter'):
    # Import APOGEE and astroNN data
    astroNN_data = import_astroNN(verbose=verbose)
    
    # Set up figure and axes
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8.5, 4.5))
    bounds = {
        'left': 0.1,
        'bottom': 0.15,
        'right': 0.85,
        'top': 0.9,
        'hspace': 0.05,
        'wspace': 0.03
    }
    plt.subplots_adjust(**bounds)
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=16, width=0.03,
                          label=r'Birth radius [kpc]', labelpad=6)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(4))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(1))
    
    for j, output in enumerate(OUTPUTS):
        if verbose: 
            print('Importing VICE multizone data from %s.vice' % output)
        vice_stars = multioutput_to_pandas(output)
        for i, absz_lim in enumerate(ABSZ_BINS):
            ax = axs[i,j]
            # Scatter plot of VICE abundances
            vice_subset = filter_multioutput_stars(vice_stars, GALR_LIM, absz_lim, 
                                                   ZONE_WIDTH)
            plot_vice_sample(ax, vice_subset, 'age', '[o/fe]', markersize=1,
                             cmap=cmap, norm=cbar.norm)
            # Median stellar ages from astroNN
            astroNN_subset = apogee_region(astroNN_data, GALR_LIM, absz_lim)
            plot_astroNN_medians(ax, astroNN_subset.copy(), ofe_lim=OFE_LIM,
                                 plot_low_count_bins=False, 
                                 low_count_cutoff=0.01, label='astroNN',
                                 markersize=4)
            # Median stellar ages from VICE
            plot_vice_medians(ax, vice_subset.copy(), ofe_lim=OFE_LIM,
                              plot_low_mass_bins=False, low_mass_cutoff=0.01,
                              label='VICE', markersize=4)
        
            # Label z-height bins
            if j == 0:
                ax.text(0.05, 0.92, r'$%s \leq |z| < %s$ kpc' % absz_lim,
                        transform=ax.transAxes, va='top', ha='left')
            
            if j == 1 and i == 0:
                ax.legend(loc='upper left', frameon=False,
                          borderaxespad=0.2, handletextpad=0.4)
        # ax.set_title(r'%s kpc $\leq |z| <$ %s kpc' % absz_lim, pad=12)
        
    # Legend
        
    # Configure axes
    for ax, label in zip(axs[0,:], LABELS):
        ax.set_title(label)
    
    for ax in axs[-1,:]:
        ax.set_xlabel('Age [Gyr]')#, labelpad=6)
    if log:
        axs[0,0].set_xscale('log')
        axs[0,0].set_xlim(AGE_LIM_LOG)
        axs[0,0].xaxis.set_major_formatter(FormatStrFormatter('%d'))
    else:
        axs[0,0].set_xlim(AGE_LIM_LINEAR)
        axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    
    for ax in axs[:,0]:
        ax.set_ylabel('[O/Fe]', labelpad=6)
    axs[0,0].set_ylim(OFE_LIM)
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    plt.savefig(paths.figures / 'presentation_age_ofe_conroy22.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_age_ofe_conroy22.py',
        description='Generate a plot of [O/Fe] vs age from VICE outputs' + \
            ' with the Conroy+ 2022 SFE for my talk at AAS 241.'
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
        default='winter',
        help='Name of colormap for color-coding VICE output (default: winter)')
    parser.add_argument('-l', '--log', action='store_true',
                        help='Plot age on a log scale')
    args = parser.parse_args()
    main(**vars(args))
