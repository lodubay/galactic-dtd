"""
Compare [O/Fe]-[Fe/H] plots for the two-infall SFH model in the inner and
outer Galaxy. Assumes a single form for the DTD.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.colors import BoundaryNorm
from matplotlib.legend_handler import HandlerTuple
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from apogee_tools import import_apogee
from ofe_feh_dtd import apogee_contours
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH, MAX_SF_RADIUS, ABSZ_BINS
import paths

FEH_LIM = (-1.4, 0.7)
OFE_LIM = (-0.15, 0.7)
GALR_BINS = [(3, 5), (11, 13)]

SFH_MODEL = 'twoinfall'
DTD_MODEL = 'plateau_width10'

def main(cmap_name='winter_r', style='paper'):
    # Set up the figure
    plt.style.use(paths.styles / f'{style}.mplstyle')
    width = ONE_COLUMN_WIDTH
    fig, axs = plt.subplots(3, 2, sharex=True, sharey=True,
                            figsize=(width, (3.5/2.5)*width))
    plt.subplots_adjust(top=0.86, right=0.92, left=0.12, bottom=0.08, 
                        wspace=0., hspace=0.)
    # Add colorbar
    birth_galr_bounds = [0, 2, 4, 6, 8, 10, 12, 14, 15.5]
    cbar = setup_colorbar(fig, cmap=cmap_name, bounds=birth_galr_bounds,
                          label=r'Birth $R_{\rm{gal}}$ [kpc]',
                          width=0.04, pad=0.02, labelpad=2)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    ism_track_color = 'k'
    ism_track_width = 0.5
    
    output_name = '/'.join(['gaussian', SFH_MODEL, DTD_MODEL, 'diskmodel'])
    # Import multioutput stars data
    mzs = MultizoneStars.from_output(output_name)
    mzs.model_uncertainty(inplace=True)
    
    for j, galr_lim in enumerate(GALR_BINS):
        axs[0,j].set_title(r'$%s\leq R_{\rm gal}<%s$ kpc' % galr_lim)
        for i in range(len(ABSZ_BINS) - 1):
            absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
            vice_subset = mzs.region(galr_lim, absz_lim)
            # Plot sample of star particle abundances
            vice_subset.scatter_plot(axs[i,j], '[fe/h]', '[o/fe]', 
                                      color='galr_origin', markersize=0.1,
                                      cmap=cmap_name, norm=cbar.norm)
            # Plot abundance tracks
            zone = int(0.5 * (galr_lim[0] + galr_lim[1]) / ZONE_WIDTH)
            zone_path = str(mzs.fullpath / ('zone%d' % zone))
            hist = vice.history(zone_path)
            axs[i,j].plot(hist['[fe/h]'], hist['[o/fe]'], c=ism_track_color, 
                          ls='-', linewidth=ism_track_width)
            # Plot APOGEE contours
            apogee_contours(axs[i,j], apogee_data, galr_lim, absz_lim)
    
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
    for i, ax in enumerate(axs[:,0]):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        ax.text(0.07, 0.93, r'$%s\leq |z| < %s$ kpc' % absz_lim, 
                va='top', ha='left', transform=ax.transAxes,
                bbox={
                    'facecolor': 'w',
                    'edgecolor': 'none',
                    'boxstyle': 'round',
                    'pad': 0.15,
                    'alpha': 1.,
                })
    # Custom legend
    custom_lines = [Line2D([0], [0], color=ism_track_color, linestyle='-', 
                           linewidth=ism_track_width),
                    (Line2D([0], [0], color='r', linestyle='-', linewidth=0.5),
                    Line2D([0], [0], color='r', linestyle='--', linewidth=0.5))]
    # legend_labels = ['Gas abundance', 'APOGEE 30% cont.', 'APOGEE 80% cont.']
    legend_labels = ['Gas abundance', 'APOGEE contours']
    axs[0,1].legend(custom_lines, legend_labels, frameon=False, edgecolor='w',
                    loc='upper right', handlelength=1.2, handletextpad=0.4,
                    borderpad=0.1, labelspacing=0.4, borderaxespad=0.5,
                    framealpha=1., 
                    handler_map={tuple: HandlerTuple(ndivide=None)})
    # Figure title with SFH and DTD models used
    fig.suptitle('Two-infall SFH + \nPlateau ($W=1$ Gyr) DTD')
    
    plt.savefig(paths.figures / 'ofe_feh_twoinfall')
    plt.close()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh_twoinfall.py',
        description='Compare [O/Fe]-[Fe/H] plots for the two-infall SFH model'+
        'in the inner and outer Galaxy. Assumes a single form for the DTD.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
