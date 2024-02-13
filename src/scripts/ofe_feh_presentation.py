"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different delay time distributions. For use in a presentation.
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from apogee_tools import import_apogee, gen_kde
from ofe_feh_dtd import apogee_contours
from utils import scatter_hist
from colormaps import paultol
from _globals import ZONE_WIDTH, TWO_COLUMN_WIDTH, MAX_SF_RADIUS, ABSZ_BINS
import paths

FEH_LIM = (-1.2, 0.6)
OFE_LIM = (-0.1, 0.65)
GALR_LIM = (7, 9)

SFH_MODEL = 'insideout'
DTD_LIST = ['powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10']
DTD_LABELS = ['Observational', 
              'Theoretical SD',
              'Theoretical DD']

CMAP_NAME = 'gray'

def main():
    # Set up the figure
    plt.style.use(paths.styles / 'presentation.mplstyle')
    fig, axs = plt.subplots(3, len(DTD_LIST), sharex=True, sharey=True,
                            figsize=(8, 6))
    plt.subplots_adjust(top=0.93, right=0.93, left=0.1, bottom=0.11, 
                        wspace=0.05, hspace=0.)
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=10, vmax=300,
                          label='# model stellar populations', labelpad=6,
                          lognorm=True)
    cbar.ax.yaxis.set_major_formatter(ScalarFormatter())
    cbar.ax.set_yticks([10, 30, 100, 300])
    # cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    # cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    # ism_track_color = 'k'
    # ism_track_width = 0.5
    
    contour_color = paultol.vibrant.colors[3]
    
    for j, dtd in enumerate(DTD_LIST):
        output_name = '/'.join(['gaussian', SFH_MODEL, dtd, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(inplace=True)
        for i in range(len(ABSZ_BINS) - 1):
            absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
            vice_subset = mzs.region(GALR_LIM, absz_lim)
            vice_sample = vice_subset.sample(10000)
            # 2D histogram + scatter plot of abundances
            # vice_subset.scatter_plot(axs[i,j], '[fe/h]', '[o/fe]', 
            #                           color='k', markersize=0.1,)
                                      # cmap=CMAP_NAME, norm=cbar.norm)
            scatter_hist(axs[i,j], vice_sample['[fe/h]'], vice_sample['[o/fe]'],
                         xlim=FEH_LIM, ylim=OFE_LIM, nbins=50, scatter_size=1,
                         log_norm=True, cmap=CMAP_NAME,
                         vmin=cbar.norm.vmin, vmax=cbar.norm.vmax)
            # print(np.nanmax(h))
            # Plot abundance tracks
            # zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
            # zone_path = str(mzs.fullpath / ('zone%d' % zone))
            # hist = vice.history(zone_path)
            # axs[i,j].plot(hist['[fe/h]'], hist['[o/fe]'], c=ism_track_color, 
            #               ls='-', linewidth=ism_track_width)
            # Plot APOGEE contours
            apogee_contours(axs[i,j], apogee_data, GALR_LIM, absz_lim,
                            linewidths=1, colors=contour_color)
    
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
    # for i, ax in enumerate(axs[:,-1]):
    #     absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
    #     ax.yaxis.set_label_position('right')
    #     ax.set_ylabel(r'$|z| = %s - %s$' % absz_lim, size=16)
    for i, ax in enumerate(axs[:,1]):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        ax.text(0.5, 0.93, r'$%s\leq |z| < %s$ kpc' % absz_lim, size=14,
                va='top', ha='center', transform=ax.transAxes,)
                # bbox={
                #     'facecolor': 'w',
                #     'edgecolor': 'none',
                #     'boxstyle': 'round',
                #     'pad': 0.15,
                #     'alpha': 1.,
                # })
    for j, ax in enumerate(axs[0]):
        ax.set_title(DTD_LABELS[j])
    # Custom legend
    axs[0,-1].text(0.95, 0.93, 'APOGEE data', color=contour_color, 
                   va='top', ha='right', transform=ax.transAxes,)
                   # bbox={
                   #     'facecolor': 'w',
                   #     'edgecolor': 'none',
                   #     'boxstyle': 'round',
                   #     'pad': 0.15,
                   #     'alpha': 1.,
                   # })
    # custom_lines = [Line2D([0], [0], color='r', linestyle='-', linewidth=1),
    #                 Line2D([0], [0], color='r', linestyle='--', linewidth=1)]
    # legend_labels = ['APOGEE 30% cont.', 'APOGEE 80% cont.']
    # axs[2, 0].legend(custom_lines, legend_labels, frameon=False, 
    #                  loc='upper left', handlelength=1, handletextpad=0.5, 
    #                  fontsize=14)
    
    plt.savefig(paths.figures / 'ofe_feh_presentation')
    plt.close()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh_presentation.py',
        description='Presentation version of the [O/Fe]-[Fe/H] DTD comparison plot',
        )
    args = parser.parse_args()
    main(**vars(args))
