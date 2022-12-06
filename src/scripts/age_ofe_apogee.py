"""
Plot a grid of [O/Fe] vs age panels for APOGEE data. Takes a minute to run.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from _globals import GALR_BINS, ABSZ_BINS
from ofe_feh_vice import setup_colorbar
from ofe_feh_apogee import read_kde, save_kde, kde2D, kde_path
import paths
from utils import import_astroNN, apogee_region, select_giants, scatter_hist, \
    get_bin_centers

NBINS = 50
GALR_BINS = GALR_BINS[:-1]
AGE_LIM = (-1, 14)
OFE_LIM = (-0.15, 0.55)

def main(verbose=True, cmap='RdPu'):
    if verbose:
        print('Importing astroNN data...')
    data = select_giants(import_astroNN())
    data.dropna(subset=['O_FE', 'ASTRONN_AGE'], inplace=True)
    if verbose:
        print('Plotting 2D histograms...')
    fig, axs = plot_scatter_hist_grid(data)
    if verbose:
        print('Plotting medians...')
    plot_contours(axs, data, cmap=cmap)
    plot_medians(axs, data)
    output = paths.figures / 'age_ofe_apogee.pdf'
    plt.savefig(output, dpi=300)
    plt.close()
    if verbose:
        print('Saved to %s.' % str(output))


def plot_scatter_hist_grid(data):
    """
    Setup a grid of panels and plot a hybrid scatterplot/2D histogram
    of APOGEE data.
    """
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1,
                          xlim=AGE_LIM, ylim=OFE_LIM)
    norm = normalize_colorbar(data, vmin=5)
    setup_colorbar(fig, 'gray', norm, label='Number of stars')
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = apogee_region(data, galr_lim, absz_lim)
            # Scatter plot / 2D histogram
            scatter_hist(ax, subset['ASTRONN_AGE'], subset['O_FE'],
                         xlim=AGE_LIM, ylim=OFE_LIM, nbins=NBINS,
                         vmin=norm.vmin, vmax=norm.vmax, cmin=norm.vmin)
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('Age [Gyr]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.9, r'$%s\leq |z| < %s$' % absz_lim,
                        transform=ax.transAxes, size=8, ha='left', va='top')
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
    return fig, axs


def plot_medians(axs, data, ofe_lim=OFE_LIM, ofe_bin_width=0.05):
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    data['OFE_BIN'] = pd.cut(data['O_FE'], ofe_bins, 
                             labels=get_bin_centers(ofe_bins))
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = apogee_region(data, galr_lim, absz_lim)
            age_grouped = subset.groupby('OFE_BIN')['ASTRONN_AGE']
            age_median = age_grouped.median()
            ax.errorbar(age_median, age_median.index, 
                        xerr=(age_median - age_grouped.quantile(0.16),
                              age_grouped.quantile(0.84) - age_median),
                        yerr=ofe_bin_width/2,
                        color='r', linestyle='none', capsize=1, elinewidth=0.5,
                        capthick=0.5, marker='o', markersize=2,
            )


def plot_contours(axs, data, cmap='Greys', linewidths=0.5):
    """
    Add contours of APOGEE abundances to axes.

    Parameters
    ----------
    axs : list of axes
    filename : str
        Path to APOGEE data file
    cmap_name : str
        Name of colormap to apply to contours
    """
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            path = kde_path(galr_lim, absz_lim, 
                            savedir='../data/APOGEE/kde/age_ofe')
            if path.exists():
                xx, yy, logz = read_kde(path)
            else:
                subset = apogee_region(data, galr_lim, absz_lim)
                xx, yy, logz = kde2D(subset['ASTRONN_AGE'], subset['O_FE'], 0.03,
                                     xbins=20j)
                # Scale by total number of stars in region
                logz += np.log(subset.shape[0])
                save_kde(xx, yy, logz, path)
            # Contour levels in log-likelihood space
            levels = np.arange(8, 11, 0.5)
            ax.contour(xx, yy, logz, levels, cmap=cmap, linewidths=linewidths)


def normalize_colorbar(data, vmin=10):
    subset = apogee_region(data,
                           galr_lim=(GALR_BINS[3], GALR_BINS[4]),
                           absz_lim=(ABSZ_BINS[0], ABSZ_BINS[1]))
    H, xedges, yedges = np.histogram2d(subset['ASTRONN_AGE'], subset['O_FE'],
                                       bins=NBINS,
                                       range=[AGE_LIM, OFE_LIM])
    norm = LogNorm(vmin=vmin, vmax=H.max())
    return norm


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
    main()
