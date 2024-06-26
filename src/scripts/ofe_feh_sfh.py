"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different star formation histories.
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.colors import BoundaryNorm
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH, MAX_SF_RADIUS
import paths
from ofe_feh_dtd import apogee_contours
from apogee_tools import import_apogee

FEH_LIM = (-1.3, 0.6)
OFE_LIM = (-0.15, 0.55)
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 0.5)

SFH_LABELS = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD_MODEL = 'exponential_timescale15'
DTD_LABEL = 'Exponential ($\\tau=1.5$ Gyr) DTD'

CMAP_NAME = 'winter_r'

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True,
                            figsize=(ONE_COLUMN_WIDTH, ONE_COLUMN_WIDTH))
    plt.subplots_adjust(top=0.88, right=0.92, wspace=0., hspace=0., left=0.11,
                        bottom=0.1)
    birth_galr_bounds = [2, 4, 6, 8, 10, 12, 14, 15.5]
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, bounds=birth_galr_bounds, 
                          label=r'Birth $R_{\rm{gal}}$ [kpc]', 
                          pad=0.02, labelpad=2, width=0.04)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    for label, ax in zip(SFH_LABELS, axs.flatten()):
        sfh = ''.join(label.split('-')).lower()
        output_name = '/'.join(['gaussian', sfh, DTD_MODEL, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(apogee_data=apogee_data, inplace=True)
        mzs.region(GALR_LIM, ABSZ_LIM, inplace=True)
        # Plot sample of star particle abundances
        mzs.scatter_plot(ax, '[fe/h]', '[o/fe]', color='galr_origin',
                          cmap=CMAP_NAME, norm=cbar.norm, markersize=0.1)
        apogee_contours(ax, apogee_data, GALR_LIM, ABSZ_LIM)
        # Plot abundance tracks
        zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
        # Import post-processed output for the given annulus
        zone_path = str(mzs.fullpath / ('zone%d' % zone))
        hist = vice.history(zone_path)
        ax.plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                linewidth=0.5)
        # Label axis
        ax.set_title(label, loc='right', x=0.93, y=0.93, va='top', pad=0)
    
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
    for ax in axs[:,0]:
        ax.set_ylabel('[O/Fe]')
    # Custom legend
    custom_lines = [Line2D([0], [0], color='k', linestyle='-', linewidth=0.5),
                    (Line2D([0], [0], color='r', linestyle='-', linewidth=0.5),
                    Line2D([0], [0], color='r', linestyle='--', linewidth=0.5))]
    # Combine contour legend entries
    legend_labels = ['Gas abundance', 'APOGEE 30/80% contours']
    axs[0,0].legend(custom_lines, legend_labels, frameon=False,
                    loc='lower left', bbox_to_anchor=(0, 1, 2, 0.1),
                    ncols=3, borderaxespad=0., handlelength=1.5, 
                    columnspacing=1.0,
                    handler_map={tuple: HandlerTuple(ndivide=None)})
    # Figure title with DTD model used
    fig.suptitle(DTD_LABEL)
    
    plt.savefig(paths.figures / 'ofe_feh_sfh', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh_sfh.py',
        description='Compare [O/Fe]-[Fe/H] plots for the Solar annulus for ' +
        'VICE outputs with different star formation histories.',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
