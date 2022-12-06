"""
Generate a grid of age vs [O/Fe] plots for multiple Galactic regions.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    sample_dataframe, get_bin_centers
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH

NBINS = 50
GALR_BINS = GALR_BINS[:-1]
AGE_LIM = (0, 14)
OFE_LIM = (-0.15, 0.55)

def main(output_name, data_dir='../data/migration', cmap='winter'):
    """
    Parameters
    ----------
    output_name : string
        Path to the
    """
    # Import multioutput stars data
    stars = multioutput_to_pandas(output_name, data_dir)
    fig, axs = plot_age_ofe_stars(stars, cmap)
    plot_post_process_tracks(output_name, axs, data_dir=data_dir)
    plot_medians(axs, stars)
    plt.savefig(paths.figures / 'age_ofe_vice.png', dpi=300)
    plt.close()


def plot_age_ofe_stars(stars, cmap):
    """
    Plot just the VICE multizone stars data
    """
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1,
                               xlim=AGE_LIM, ylim=OFE_LIM)
    norm = normalize_colorbar(stars)
    cax = setup_colorbar(fig, cmap, norm, label=r'Birth $R_{\rm{Gal}}$ [kpc]')
    cax.yaxis.set_minor_locator(MultipleLocator(0.5))
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
                                              ZONE_WIDTH)
            sample = sample_dataframe(subset, 10000)
            # Scatter plot of random sample of stellar particles
            ax.scatter(sample['age'], sample['[o/fe]'], s=0.1,
                       c=sample['zone_origin'] * ZONE_WIDTH, cmap=cmap,
                       norm=norm, rasterized=True, edgecolor='none')
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('Age [Gyr]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.85, r'$%s\leq |z| < %s$' % absz_lim,
                        transform=ax.transAxes, size=8)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
    return fig, axs


def plot_post_process_tracks(output_name, axs, data_dir='../data/migration'):
    """
    Plot post-process abundance track for the mean annulus in each panel.

    Parameters
    ----------
    output_name : str
        Path to the VICE output directory
    axs : list of axes
    data_dir : str, optional
        Parent directory for all migration outputs
    """
    # Convert annulus from kpc to zone number
    nzones = len(GALR_BINS)-1
    mean_galr = [(GALR_BINS[i] + GALR_BINS[i+1]) / 2 for i in range(nzones)]
    zones = [int(mean_galr[i] / ZONE_WIDTH) for i in range(nzones)]
    for i, zone in enumerate(zones):
        # Import post-processed output for the given annulus
        # post_process_path = Path(data_dir) / Path(output_name + '.vice') / ('zone%d' % zone)
        post_process_path = Path(data_dir) / Path(output_name.replace(
            'diffusion', 'post-process') + '.vice') / Path('zone%d' % zone)
        post_process_hist = vice.history(str(post_process_path))
        for ax in axs[:,i].flatten():
            ax.plot(post_process_hist['lookback'], post_process_hist['[o/fe]'],
                    c='k', ls='-', linewidth=1)


def plot_medians(axs, stars, ofe_lim=OFE_LIM, ofe_bin_width=0.05):
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    stars['odf_bin'] = pd.cut(stars['[o/fe]'], ofe_bins, 
                             labels=get_bin_centers(ofe_bins))
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
                                              ZONE_WIDTH)
            age_grouped = subset.groupby('odf_bin')['age']
            age_median = age_grouped.median()
            ax.errorbar(age_median, age_median.index, 
                        xerr=(age_median - age_grouped.quantile(0.16),
                              age_grouped.quantile(0.84) - age_median),
                        yerr=ofe_bin_width/2,
                        color='k', linestyle='none', capsize=1, elinewidth=0.5,
                        capthick=0.5, marker='o', markersize=2,
            )


def normalize_colorbar(stars):
    """
    Calculate the colorbar normalization for the full migration output data
    before subsampling.

    Parameters
    ----------
    stars : pandas DataFrame
        VICE migration output

    Returns
    -------
    norm : instance of matplotlib.colors.Normalize
    """
    # Filter to full range to be plotted
    subset = filter_multioutput_stars(stars,
                                      galr_lim=(GALR_BINS[0], GALR_BINS[-1]),
                                      absz_lim=(ABSZ_BINS[0], ABSZ_BINS[-1]),
                                      zone_width=ZONE_WIDTH)
    norm = Normalize(vmin=0, vmax=subset['zone_origin'].max() * ZONE_WIDTH)
    return norm


def setup_colorbar(fig, cmap, norm, label=''):
    """
    Configure colorbar given a colormap and a normalization.

    Parameters
    ----------
    fig
    cmap
    norm

    Returns
    -------
    cax : colorbar axis
    """
    # Colorbar axis
    plt.subplots_adjust(right=0.92, left=0.06, bottom=0.07, top=0.95,
                        wspace=0.05, hspace=0.05)
    cax = plt.axes([0.93, 0.07, 0.02, 0.88])
    # Add colorbar
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label(label)
    return cax


def setup_axes(rows, cols, width=8, xlim=None, ylim=None):
    """
    Set up a blank grid of axes plus a colorbar axis.

    Parameters
    ----------
    rows : int
        Number of rows of axes
    cols : int
        Number of columns of axes
    width : float, optional
        Width of the figure in inches. The default is 8 in.
    xlim : tuple or None, optional
        Limits of x-axis for all axes
    ylim : tuple or None, optional
        Limits of y-axis for all axes

    Returns
    -------
    fig : matplotlib figure
    axs : list of axes
    cax : axis object for colorbar
    """
    fig, axs = plt.subplots(rows, cols, figsize=(width, (width/cols)*rows),
                            sharex=True, sharey=True)
    # Configure axis limits and ticks (will be applied to all axes)
    axs[0,0].set_xlim(xlim)
    axs[0,0].set_ylim(ylim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.1))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.02))
    return fig, axs


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))
