"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different delay time distributions. For use in a presentation.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import ScalarFormatter
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from apogee_tools import import_apogee
from ofe_feh_dtd import apogee_contours
from utils import scatter_hist
from colormaps import paultol
# from _globals import ABSZ_BINS
import paths

FEH_LIM = (-1.2, 0.6)
OFE_LIM = (-0.1, 0.65)
GALR_LIM = (7, 9)
ABSZ_BINS = [(0, 0.5), (1, 2)]

SFH_MODEL = 'insideout'
DTD_LIST = ['powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10']
DTD_LABELS = ['Observational', 
              'Theoretical\nWD+star',
              'Theoretical\nWD+WD']
ROW_LABELS = ['Far from midplane', 'Near midplane']

CMAP_NAME = 'gray'

def main(contours=False):
    # Set up the figure
    plt.style.use(paths.styles / 'presentation.mplstyle')
    fig, axs = plt.subplots(len(ABSZ_BINS), len(DTD_LIST), sharex=True, 
                            sharey=True, figsize=(8, 5))
    plt.subplots_adjust(top=0.87, right=0.93, left=0.1, bottom=0.13, 
                        wspace=0.05, hspace=0.)
    # Vertical colorbar for 2d density histogram
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=10, vmax=300,
                          label='# model stellar populations', labelpad=6,
                          lognorm=True)
    cbar.ax.yaxis.set_major_formatter(ScalarFormatter())
    cbar.ax.set_yticks([10, 30, 100, 300])
    
    if contours:
        print('Importing APOGEE data...')
        apogee_data = import_apogee()
    
    # contour_color = paultol.vibrant.colors[3]
    contour_color = 'r'
    
    print('Plotting...')
    for j, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH_MODEL, dtd, 'diskmodel'])
        print('\t%s [%s/%s]' % (output_name, j+1, len(DTD_LIST)))
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(inplace=True)
        for i in range(len(ABSZ_BINS)):
            # absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
            absz_lim = ABSZ_BINS[-(i+1)]
            vice_subset = mzs.region(GALR_LIM, absz_lim)
            vice_sample = vice_subset.sample(10000)
            # 2D histogram + scatter plot of abundances
            scatter_hist(axs[i,j], vice_sample['[fe/h]'], vice_sample['[o/fe]'],
                         xlim=FEH_LIM, ylim=OFE_LIM, nbins=50, scatter_size=1,
                         log_norm=True, cmap=CMAP_NAME,
                         vmin=cbar.norm.vmin, vmax=cbar.norm.vmax)
            # Plot APOGEE contours
            if contours:
                apogee_contours(axs[i,j], apogee_data, GALR_LIM, absz_lim,
                                linewidths=1.5, colors=contour_color)
    
    # Set x-axis ticks
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # Set y-axis ticks
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    # Set axis limits
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    # Set axis labels
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    for i, ax in enumerate(axs[:,0]):
        ax.set_ylabel('[O/Fe]', labelpad=6)
    # Label rows with z-height bounds
    for i, ax in enumerate(axs[:,0]):
        # absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        # ax.text(0.5, 0.93, r'$%s\leq |z| < %s$ kpc' % absz_lim, size=14,
        #         va='top', ha='center', transform=ax.transAxes)
        ax.text(0.5, 0.93, ROW_LABELS[i],
                va='top', ha='center', transform=ax.transAxes)
    # Axis column titles
    for j, ax in enumerate(axs[0]):
        ax.set_title(DTD_LABELS[j])
    # Label APOGEE data (no legend)
    if contours:
        axs[0,-1].text(0.95, 0.93, 'Observed stars', color=contour_color, 
                       va='top', ha='right', transform=ax.transAxes)
    
    fname = 'presentation/ofe_feh_dtd'
    if not contours:
        fname += '_nocontours'
    plt.savefig(paths.figures / fname)
    plt.close()
    print('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_ofe_feh.py',
        description='Presentation version of the [O/Fe]-[Fe/H] DTD comparison plot',
        )
    parser.add_argument('-c', '--contours', action='store_true',
                        help='Overplot APOGEE data contours in each panel.')
    args = parser.parse_args()
    main(**vars(args))
