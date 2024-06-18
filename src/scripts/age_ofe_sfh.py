"""
Compare [O/Fe]-Age plots for the Solar annulus for VICE outputs with
different star formation histories.
"""

import argparse
import numpy as np
from apogee_tools import import_apogee, apogee_region
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from multizone_stars import MultizoneStars
from age_ofe import plot_vice_medians, plot_astroNN_medians
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH, MAX_SF_RADIUS
import paths

SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
DTD_MODEL = 'exponential_timescale15'
DTD_LABEL = 'Exponential ($\\tau=1.5$ Gyr) DTD'
LABEL_LIST = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
AGE_COL = 'LATENT_AGE'
AGE_LABEL = 'L23'#'Leung et al.\n(2023)'
AGE_LIM = (0.3, 20)
OFE_LIM = (-0.15, 0.5)
CMAP_NAME = 'viridis'
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 0.5)

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    apogee_data = import_apogee()
    apogee_subset = apogee_region(apogee_data, GALR_LIM, ABSZ_LIM)
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True,
                            figsize=(ONE_COLUMN_WIDTH, 0.92*ONE_COLUMN_WIDTH))
    plt.subplots_adjust(top=0.92, right=0.85, wspace=0., hspace=0., bottom=0.11)
    # Define colorbar axis
    height = fig.subplotpars.top - fig.subplotpars.bottom - 0.06
    cax = plt.axes([fig.subplotpars.right + 0.02, fig.subplotpars.bottom, 
                    0.04, height])
    # Add colorbar
    # norm = Normalize(vmin=-1.3, vmax=0.3)
    bounds = np.arange(-1.3, 0.5, 0.2)
    cmap = plt.get_cmap(CMAP_NAME)
    norm = BoundaryNorm(bounds, cmap.N, extend='both')
    cbar = fig.colorbar(ScalarMappable(norm, CMAP_NAME), cax)
    # align title to colorbar bounding box
    bbox = cbar.ax.get_window_extent()
    x, _ = cbar.ax.transAxes.inverted().transform([bbox.x0, bbox.y0])
    cbar.ax.set_title('[Fe/H]', ha='left', x=x)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(0.5))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    
    for ax, sfh, label in zip(axs.flatten(), SFH_LIST, LABEL_LIST):
        output_name = '/'.join(['gaussian', sfh, DTD_MODEL, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data=apogee_data, inplace=True)
        mzs.region(GALR_LIM, ABSZ_LIM, inplace=True)
        # Plot sample of star particle abundances
        mzs.scatter_plot(ax, 'age', '[o/fe]', color='[fe/h]',
                          cmap=CMAP_NAME, norm=cbar.norm)
        plot_astroNN_medians(ax, apogee_subset, age_col=AGE_COL, 
                              label=AGE_LABEL, 
                              plot_low_count_bins=False)
        plot_vice_medians(ax, mzs.stars, label='Model',
                          plot_low_mass_bins=False)
        # Label axis
        ax.set_title(label, loc='left', x=0.07, y=0.93, va='top', pad=0)
    
    # Set x-axis scale and ticks
    axs[0,0].set_xlim(AGE_LIM)
    axs[0,0].set_xscale('log')
    axs[0,0].xaxis.set_major_formatter(FormatStrFormatter('%d'))
    
    # Set y-axis ticks
    axs[0,0].set_ylim(OFE_LIM)
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    # Axis labels
    for ax in axs[-1]:
        ax.set_xlabel('Age [Gyr]')
    for i, ax in enumerate(axs[:,0]):
        ax.set_ylabel('[O/Fe]', labelpad=2)
        
    # Legend
    axs[0,0].legend(loc='upper left', frameon=False, 
                    bbox_to_anchor=(0.02, 0.89), handlelength=0.7)
    
    # Figure title
    fig.suptitle(DTD_LABEL)
    
    fig.savefig(paths.figures / 'age_ofe_sfh', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='age_feh_sfh.py',
        description='Compare Age-[O/Fe] plots for the Solar annulus for ' +
        'VICE outputs with different star formation histories.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
