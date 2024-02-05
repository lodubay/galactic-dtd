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
from _globals import ZONE_WIDTH, TWO_COLUMN_WIDTH, MAX_SF_RADIUS, ABSZ_BINS
import paths

FEH_LIM = (-1.3, 0.6)
OFE_LIM = (-0.1, 0.6)
GALR_LIM = (11, 13)

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
        for i in range(len(ABSZ_BINS) - 1):
            absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
            vice_subset = mzs.region(GALR_LIM, absz_lim)
            # Plot sample of star particle abundances
            vice_subset.scatter_plot(axs[i,j], '[fe/h]', '[o/fe]', 
                                      color='galr_origin', markersize=0.1,
                                      cmap=CMAP_NAME, norm=cbar.norm)
            # Plot abundance tracks
            zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
            zone_path = str(mzs.fullpath / ('zone%d' % zone))
            hist = vice.history(zone_path)
            axs[i,j].plot(hist['[fe/h]'], hist['[o/fe]'], c=ism_track_color, 
                          ls='-', linewidth=ism_track_width)
            # Plot APOGEE contours
            apogee_contours(axs[i,j], apogee_data, GALR_LIM, absz_lim)
    
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
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        ax.text(0.5, 0.93, r'$%s\leq |z| < %s$ kpc' % absz_lim, 
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
    axs[2, 0].legend(custom_lines, legend_labels, frameon=False, 
                     loc='upper left', handlelength=0.6, handletextpad=0.4)
    
    plt.savefig(paths.figures / 'ofe_feh_dtd_outer')
    plt.close()


def apogee_contours(ax, apogee_data, galr_lim=(0, 20), absz_lim=(0, 3)):
    xx, yy, logz = gen_kde(apogee_data, bandwidth=0.03,
                           galr_lim=GALR_LIM, absz_lim=absz_lim)
    # scale the linear density to the max value
    scaled_density = np.exp(logz) / np.max(np.exp(logz))
    # contour levels at 1 and 2 sigma
    levels = get_levels(scaled_density)
    # levels = np.exp(-0.5 * np.array([2, 1])**2)
    cs = ax.contour(xx, yy, scaled_density, levels, colors='r',
                    linewidths=0.5, linestyles=['--', '-'])


def get_levels(scaled_density, enclosed=[0.8, 0.3]):
    """
    Calculate the contour levels which contain the given enclosed probabilities.
    """
    levels = []
    l = 0.
    i = 0
    while l < 1 and i < len(enclosed):
        frac_enclosed = np.sum(scaled_density[scaled_density > l]) / np.sum(scaled_density)
        if frac_enclosed <= enclosed[i] + 0.01:
            levels.append(l)
            i += 1
        l += 0.01
    return levels
    


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
