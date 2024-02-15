"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different delay time distributions.
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from apogee_tools import import_apogee, gen_kde
from ofe_feh_dtd import apogee_contours
from _globals import ZONE_WIDTH, TWO_COLUMN_WIDTH, MAX_SF_RADIUS
import paths

FEH_LIM = (-1.3, 0.6)
OFE_LIM = (-0.1, 0.6)
ABSZ_LIM = (0.5, 1)
GALR_BINS = [(11, 13), (7, 9), (3, 5)]

SFH_MODEL = 'insideout'
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
DTD_LABELS = ['Two-population', 
              'Power-law\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1$ Gyr)',
              'Triple-system']

CMAP_NAME = 'winter'

def main(style='paper'):
    # Set up the figure
    plt.style.use(paths.styles / f'{style}.mplstyle')
    width = TWO_COLUMN_WIDTH
    fig, axs = plt.subplots(3, 5, sharex=True, sharey=True,
                            figsize=(width, 3/5*width))
    plt.subplots_adjust(top=0.91, right=0.98, left=0.06, bottom=0.08, 
                        wspace=0., hspace=0.)
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=0, vmax=MAX_SF_RADIUS,
                          label=r'Birth $R_{\rm{gal}}$ [kpc]')
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    ism_track_color = 'k'
    ism_track_width = 0.5
    
    for j, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH_MODEL, dtd, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(inplace=True)
        for i, galr_lim in enumerate(GALR_BINS):
            vice_subset = mzs.region(galr_lim, ABSZ_LIM)
            # Plot sample of star particle abundances
            vice_subset.scatter_plot(axs[i,j], '[fe/h]', '[o/fe]', 
                                      color='galr_origin', markersize=0.1,
                                      cmap=CMAP_NAME, norm=cbar.norm)
            # Plot abundance tracks
            zone = int(0.5 * (galr_lim[0] + galr_lim[1]) / ZONE_WIDTH)
            zone_path = str(mzs.fullpath / ('zone%d' % zone))
            hist = vice.history(zone_path)
            axs[i,j].plot(hist['[fe/h]'], hist['[o/fe]'], c=ism_track_color, 
                          ls='-', linewidth=ism_track_width)
            # Plot APOGEE contours
            apogee_contours(axs[i,j], apogee_data, galr_lim, ABSZ_LIM)
    
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
        ax.set_ylabel('[O/Fe]', labelpad=2)
    for i, ax in enumerate(axs[:,2]):
        galr_lim = GALR_BINS[i]
        ax.text(0.5, 0.93, r'$%s\leq R_{\rm gal} < %s$ kpc' % galr_lim, 
                va='top', ha='center', transform=ax.transAxes,
                bbox={
                    'facecolor': 'w',
                    'edgecolor': 'none',
                    'boxstyle': 'round',
                    'pad': 0.15,
                    'alpha': 1.,
                })
    for j, ax in enumerate(axs[0]):
        ax.set_title(DTD_LABELS[j])
    # Custom legend
    custom_lines = [Line2D([0], [0], color=ism_track_color, linestyle='-', 
                           linewidth=ism_track_width),
                    Line2D([0], [0], color='r', linestyle='-', linewidth=0.5),
                    Line2D([0], [0], color='r', linestyle='--', linewidth=0.5)]
    legend_labels = ['Gas abundance', 'APOGEE 30% cont.', 'APOGEE 80% cont.']
    axs[0, 0].legend(custom_lines, legend_labels, frameon=False, 
                     loc='upper left', handlelength=0.6, handletextpad=0.4)
    
    plt.savefig(paths.figures / 'ofe_feh_dtd_Rgal')
    plt.close()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh_dtd.py',
        description='Compare [O/Fe]-[Fe/H] plots for the Solar annulus for ' +
        'VICE outputs with different delay time distributions.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
