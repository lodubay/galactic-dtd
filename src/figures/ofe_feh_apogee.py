"""
Plot a grid of [O/Fe] vs [Fe/H] panels for APOGEE data
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from sklearn.neighbors import KernelDensity
from ofe_feh_vice import GALR_BINS, ABSZ_BINS, FEH_LIM, OFE_LIM
from ofe_feh_vice import setup_axes, setup_colorbar

global NBINS
NBINS = 20

def main(filename='../data/APOGEE/dr17_cut_data.csv', cmap_name='magma'):
    print('Importing data')
    data = pd.read_csv(Path(filename))
    print('Scatterplots and histograms')
    fig, axs = plot_scatter_hist_grid(data, contour_cmap=cmap_name)
    plt.savefig('ofe_feh_apogee.pdf', dpi=300)
    plt.close()


def plot_scatter_hist_grid(data, contour_cmap='magma'):
    """
    Setup a grid of panels and plot a hybrid scatterplot/2D histogram
    of APOGEE data.
    """
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1,
                          xlim=FEH_LIM, ylim=OFE_LIM)
    norm = normalize_colorbar(data)
    setup_colorbar(fig, 'gray', norm, label='Number of stars')
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            print('\tPanel %s, %s' % (i, j))
            print('\t\tSlicing data')
            subset = apogee_region(data, galr_lim, absz_lim)
            # Scatter plot / 2D histogram
            print('\t\tScatter plot')
            scatter_hist(ax, subset['FE_H'], subset['O_FE'],
                         xlim=FEH_LIM, ylim=OFE_LIM, nbins=NBINS,
                         vmin=norm.vmin, vmax=norm.vmax)
            # Contour plot
            print('\t\tContourPlot')
            x, y, z = kde2D(subset['FE_H'], subset['O_FE'], 0.05)
            ax.contour(x, y, z, cmap=contour_cmap)
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.55, 0.85, r'$%s\leq |z| < %s$' % absz_lim,
                        transform=ax.transAxes)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
    return fig, axs


def kde2D(x, y, bandwidth, xbins=50j, ybins=50j, **kwargs):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins,
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)


def normalize_colorbar(data):
    subset = apogee_region(data,
                           galr_lim=(GALR_BINS[2], GALR_BINS[3]),
                           absz_lim=(ABSZ_BINS[0], ABSZ_BINS[1]))
    H, xedges, yedges = np.histogram2d(subset['FE_H'], subset['O_FE'],
                                       bins=NBINS,
                                       range=[FEH_LIM, OFE_LIM])
    norm = LogNorm(vmin=10, vmax=H.max())
    return norm


def apogee_region(data, galr_lim=(0, 20), absz_lim=(0, 5)):
    """
    Slice APOGEE data within a given Galactic region of radius and z-height.

    Parameters
    ----------
    stars : pandas DataFrame
        Output from stars_dataframe()
    galr_lim : tuple
        Minimum and maximum Galactic radius in kpc
    absz_lim : tuple
        Minimum and maximum of the absolute value of z-height in kpc
    zone_width : float
        Width of each simulation zone in kpc

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of stellar parameters
    """
    galr_min, galr_max = galr_lim
    absz_min, absz_max = absz_lim
    # Select subset
    subset = data[(data['ASTRONN_GALR'] >= galr_min) &
                  (data['ASTRONN_GALR'] < galr_max) &
                  (data['ASTRONN_GALZ'].abs() >= absz_min) &
                  (data['ASTRONN_GALZ'].abs() < absz_max)]
    subset.reset_index(inplace=True)
    return subset.dropna(subset='O_FE')


def scatter_hist(ax, x, y, xlim=None, ylim=None, log_norm=True, cmap='gray',
                 cmin=10, vmin=None, vmax=None, nbins=50, color='k',
                 rasterized=True):
    """
    Generate a scatter plot and overlayed 2D histogram for dense data.

    Parameters
    ----------
    ax : matplotlib.axis.Axes
        Axes object on which to plot the data.
    x : array-like
        Horizontal coordinates of the data points.
    y : array-like
        Vertical coordinates of the data points.
    xlim : float, optional
        Bounds for x-axis. The default is None.
    ylim : float, optional
        Bounds for y-axis. The default is None.
    log_norm : bool, optional
        Shade the 2D histogram on a logarithmic scale. The default is True.
    cmap : str, optional
        Colormap for 2D histogram. The default is'gray'.
    cmin : int, optional
        Minimum counts per bin; any number below this will show individual points.
        The default is 10.
    vmin : float or None, optional
        Value to map to minimum of histogram normalization. The default is None.
    vmax : float or None, optional
        Value to map to maximum of histogram normalization. The default is None.
    nbins : int or tuple of ints, optional
        Number of histogram bins. If a tuple, presumed to be (xbins, ybins).
        The default is 50.
    color : str, optional
        Color of individual points. The default is 'k'.
    rasterized : bool, optional [default: True]
        Whether to rasterize the scattered points
    """
    # Set automatic plot bounds
    if not xlim:
        xlim = (np.min(x), np.max(x))
    if not ylim:
        ylim = (np.min(y), np.max(y))
    # Set bin edges
    if type(nbins) == 'tuple':
        xbins, ybins = nbins
    else:
        xbins = ybins = nbins
    xbins = np.linspace(xlim[0], xlim[1], num=xbins, endpoint=True)
    ybins = np.linspace(ylim[0], ylim[1], num=ybins, endpoint=True)
    # Histogram normalization
    if log_norm:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    # Plot
    ax.scatter(x, y, c=color, s=0.5, rasterized=rasterized)
    ax.hist2d(x, y, bins=[xbins, ybins], cmap=cmap, norm=norm, cmin=cmin)


if __name__ == '__main__':
    main()
