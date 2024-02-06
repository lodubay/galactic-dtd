"""
Compare plots of [O/Fe] vs age in a single Galactic region for various DTDs.
"""

import argparse
from apogee_tools import import_apogee, apogee_region
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from multizone_stars import MultizoneStars
from age_ofe import plot_vice_medians, plot_astroNN_medians
from scatter_plot_grid import setup_colorbar
from _globals import ZONE_WIDTH, TWO_COLUMN_WIDTH, MAX_SF_RADIUS, ABSZ_BINS
import paths

SFH_MODEL = 'earlyburst'
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
LABEL_LIST = ['Two-population', 
              'Power law\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1.0$ Gyr)', 
              'Triple-system']
AGE_SOURCE = 'L23'
AGE_COL = 'LATENT_AGE'
AGE_LABEL = 'L23'#'Leung et al.\n(2023)'
AGE_LIM = (0.3, 20)
OFE_LIM = (-0.15, 0.55)
GALR_LIM = (3, 5)
CMAP_NAME = 'viridis'

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    width = TWO_COLUMN_WIDTH
    fig, axs = plt.subplots(3, 5, sharex=True, sharey=True,
                            figsize=(width, 3/5*width))
    plt.subplots_adjust(top=0.91, right=0.98, left=0.06, bottom=0.08, 
                        wspace=0., hspace=0.)
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=-1.3, vmax=0.3)
    # align title to colorbar bounding box
    bbox = cbar.ax.get_window_extent()
    x, _ = cbar.ax.transAxes.inverted().transform([bbox.x0, bbox.y0])
    cbar.ax.set_title('[Fe/H]', ha='left', x=x)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(0.5))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    
    apogee_data = import_apogee()
    
    for j, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH_MODEL, dtd, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(inplace=True)
        for i in range(len(ABSZ_BINS) - 1):
            absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
            apogee_subset = apogee_region(apogee_data, GALR_LIM, absz_lim)
            vice_subset = mzs.region(GALR_LIM, absz_lim)
            # Plot sample of star particle abundances
            vice_subset.scatter_plot(axs[i,j], 'age', '[o/fe]', 
                                     color='[fe/h]',
                                     cmap=CMAP_NAME, norm=cbar.norm)
            plot_astroNN_medians(axs[i,j], apogee_subset, age_col=AGE_COL, 
                                 label=AGE_LABEL, plot_low_count_bins=False)
            plot_vice_medians(axs[i,j], vice_subset.stars, label='Model',
                              plot_low_mass_bins=False)
    
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
    for i, ax in enumerate(axs[:,2]):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        ax.text(0.07, 0.93, r'$%s\leq |z| < %s$ kpc' % absz_lim, 
                va='top', transform=ax.transAxes)
    for j, ax in enumerate(axs[0]):
        ax.set_title(LABEL_LIST[j])
        
    # Legend
    # axs[0,0].legend(loc='upper left', frameon=False, 
    #                 bbox_to_anchor=(0.02, 0.89), handlelength=0.7)
    axs[0,-1].legend(loc='upper left', frameon=False, handlelength=0.7)
    
    fig.savefig(paths.figures / 'age_ofe_dtd_inner')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='age_ofe_dtd.py',
        description='Compare age-[O/Fe] plots for the Solar annulus for ' +
        'VICE outputs with different delay time distributions.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
