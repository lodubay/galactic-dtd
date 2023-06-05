"""
Plot a grid of [O/Fe] vs [Fe/H] panels for APOGEE data. Takes a minute to run.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from ofe_feh_vice import GALR_BINS, ABSZ_BINS, FEH_LIM, OFE_LIM
from ofe_feh_vice import setup_axes, setup_colorbar
import paths
from apogee_tools import import_apogee, apogee_region
from utils import scatter_hist, kde2D

global NBINS
NBINS = 50
DEFAULT_LINESTYLES = ['dotted', 'dashed', 'solid']

def main(cmap='winter', verbose=True):
    if verbose:
        print('Importing allStar data...')
    data = import_apogee()
    if verbose:
        print('Plotting 2D histograms...')
    fig, axs = plot_scatter_hist_grid(data)
    if verbose:
        print('Plotting contours...')
    plot_contours_grid(axs, data, cmap=cmap, colors=None)
    output = paths.figures / 'ofe_feh_apogee.png'
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
                          xlim=FEH_LIM, ylim=OFE_LIM)
    norm = normalize_colorbar(data)
    setup_colorbar(fig, 'gray', norm, label='Number of stars')
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = apogee_region(data, galr_lim, absz_lim)
            # Scatter plot / 2D histogram
            scatter_hist(ax, subset['FE_H'], subset['O_FE'],
                         xlim=FEH_LIM, ylim=OFE_LIM, nbins=NBINS,
                         vmin=norm.vmin, vmax=norm.vmax)
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.1, r'$%s\leq |z| < %s$' % absz_lim,
                        transform=ax.transAxes, size=8, va='bottom', ha='left')
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
    return fig, axs


def plot_contours_grid(axs, data, bandwidth=0.02, cmap=None, colors='k', 
                  linewidths=0.5, linestyles=DEFAULT_LINESTYLES):
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
            plot_contours(ax, data, bandwidth=bandwidth, absz_lim=absz_lim,
                          galr_lim=galr_lim, cmap=cmap, colors=colors,
                          linewidths=linewidths, linestyles=linestyles)
            

def plot_contours(ax, data, bandwidth=0.02, absz_lim=(0, 5), galr_lim=(0, 20),
                  overwrite=False, **kwargs):
    """
    Plot APOGEE density contours in the specified region.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object on which to draw contours.
    data : pandas.DataFrame
        APOGEE data containing columns 'FE_H' and 'O_FE'.
    bandwidth : float
        Kernel density estimate bandwidth. A larger number will produce
        smoother contour lines. The default is 0.02.
    absz_lim : tuple
        Limits on absolute Galactic z-height in kpc. The default is (0, 5).
    galr_lim : tuple
        Limits on Galactocentric radius in kpc. The default is (0, 20).
    overwrite : bool
        If True, force re-generate the 2D KDE and save the output.
    **kwargs passed to matplotlib.axes.Axes.contour().
    
    Returns
    -------
    contours: matplotlib.contour.QuadContourSet
    """
    xx, yy, logz = gen_kde(data, bandwidth=bandwidth, absz_lim=absz_lim, 
                           galr_lim=galr_lim, overwrite=overwrite)
    # scale the linear density to the max value
    scaled_density = np.exp(logz) / np.max(np.exp(logz))
    # contour levels at 1, 2, and 3 sigma
    levels = np.exp(-0.5 * np.array([3, 2, 1])**2)
    contours = ax.contour(xx, yy, scaled_density, levels, **kwargs)
    return contours


def gen_kde(data, bandwidth=0.02, absz_lim=(0, 5), galr_lim=(0, 20), 
            overwrite=False):
    """
    Generate kernel density estimate (KDE) of APOGEE data, or import previously
    saved KDE if it already exists.
    
    Parameters
    ----------
    data : pandas.DataFrame
        APOGEE data containing columns 'FE_H' and 'O_FE'.
    bandwidth : float
        Kernel density estimate bandwidth. A larger number will produce
        smoother contour lines. The default is 0.02.
    absz_lim : tuple
        Limits on absolute Galactic z-height in kpc. The default is (0, 5).
    galr_lim : tuple
        Limits on Galactocentric radius in kpc. The default is (0, 20).
    overwrite : bool
        If True, force re-generate the 2D KDE and save the output.
    
    Returns
    -------
    xx, yy, logz: tuple of numpy.array
        Outputs of kde2D()
    """    
    # Path to save 2D KDE for faster plot times
    path = kde_path(galr_lim, absz_lim, savedir=paths.data/'APOGEE/kde/ofe_feh/')
    if path.exists() and not overwrite:
        xx, yy, logz = read_kde(path)
    else:
        # Limit to specified Galactic region
        subset = apogee_region(data, galr_lim, absz_lim)
        subset = subset.copy().dropna(axis=0, subset=['FE_H', 'O_FE'])
        xx, yy, logz = kde2D(subset['FE_H'], subset['O_FE'], bandwidth)
        save_kde(xx, yy, logz, path)
    return xx, yy, logz


def read_kde(path):
    """
    Read a text file generated by save_kde()
    """
    arr2d = np.genfromtxt(path)
    nrows = int(arr2d.shape[0]/3)
    xx = arr2d[:nrows]
    yy = arr2d[nrows:2*nrows]
    logz = arr2d[2*nrows:]
    return xx, yy, logz


def save_kde(xx, yy, logz, path):
    """
    Generate a text file containing the KDE of the given region along with its
    corresponding x and y coordinates.
    """
    with open(path, 'w') as f:
        for arr in [xx, yy, logz]:
            f.write('#\n')
            np.savetxt(f, arr)


def kde_path(galr_lim, absz_lim, savedir=paths.data/'APOGEE/kde'):
    """
    Generate file name for the KDE of the given region.
    """
    filename = 'r%s-%s_z%s-%s.dat' % (galr_lim[0], galr_lim[1],
                                      absz_lim[0], absz_lim[1])
    return Path(savedir) / filename


def normalize_colorbar(data):
    subset = apogee_region(data,
                           galr_lim=(GALR_BINS[2], GALR_BINS[3]),
                           absz_lim=(ABSZ_BINS[0], ABSZ_BINS[1]))
    H, xedges, yedges = np.histogram2d(subset['FE_H'], subset['O_FE'],
                                       bins=NBINS,
                                       range=[FEH_LIM, OFE_LIM])
    norm = LogNorm(vmin=10, vmax=H.max())
    return norm


if __name__ == '__main__':
    main()
