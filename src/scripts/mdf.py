"""
Plot metallicity distribution functions (MDFs) of [Fe/H] binned by radius.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None # default='warn'
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from sklearn.neighbors import KernelDensity
import vice
from _globals import ABSZ_BINS, ZONE_WIDTH
from ofe_feh_apogee import apogee_region
from utils import multioutput_to_pandas, filter_multioutput_stars, import_allStar
import paths

global FEH_LIM
global GALR_BINS

FEH_LIM = (-1.1, 0.6)
BIN_WIDTH = 0.01
SMOOTH_WIDTH = 0.2
GALR_BINS = [3, 5, 7, 9, 11, 13, 15]

def main(migration, evolution, RIa, cmap_name='plasma_r'):
    output = '%s/%s/%s' % (migration, evolution, RIa)
    stars = multioutput_to_pandas(output, paths.data / 'migration')
    data = import_allStar()

    fig, axs = setup_axes(xlim=FEH_LIM)
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_colorbar(fig, cmap, norm, label=r'Galactocentric radius [kpc]')

    rmin, rmax = GALR_BINS[0], GALR_BINS[-2]
    colors = cmap([(r-rmin)/(rmax-rmin) for r in GALR_BINS[:-1]])
    for i in range(len(ABSZ_BINS)-1):
        absz_lim = ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i]
        axs[i,0].set_ylabel(r'$|z| = %s - %s$' % tuple(absz_lim))
        for j in range(len(GALR_BINS)-1):
            galr_lim = GALR_BINS[j:j+2]
            # Plot VICE post-process in left panels
            vice_subset = filter_multioutput_stars(stars,
                                                   galr_lim, absz_lim,
                                                   ZONE_WIDTH, min_mass=0)
            vice_mdf, bins = gen_mdf(vice_subset, bin_width=BIN_WIDTH)
            vice_smooth = box_smooth(vice_mdf, bins, SMOOTH_WIDTH)
            x_plot = bins[:-1] + BIN_WIDTH/2
            axs[i,0].plot(x_plot, vice_smooth, color=colors[j], linewidth=1)

            # Plot APOGEE in right panels
            apogee_subset = apogee_region(data, galr_lim, absz_lim)
            apogee_mdf, _ = np.histogram(apogee_subset['FE_H'], bins=bins, density=True)
            apogee_smooth = box_smooth(apogee_mdf, bins, SMOOTH_WIDTH)
            axs[i,1].plot(x_plot, apogee_smooth, color=colors[j], linewidth=1)

    ylim = axs[0,0].get_ylim()
    axs[0,0].set_ylim((1.5*ylim[0], None))
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    # for ax in axs[:,0]:
        # ax.set_ylabel('PDF')
    axs[0,0].set_title('/'.join((evolution, RIa)))
    # axs[0,0].set_title('VICE')
    axs[0,1].set_title('APOGEE DR17')
    plt.savefig(paths.figures / 'feh_mdf.png', dpi=300)
    plt.close()


def gen_mdf(stars, col='[fe/h]', range=None, bin_width=0.05):
    """
    Generate a metallicity distribution function (MDF) from a VICE multizone run.

    Parameters
    ----------
    stars : DataFrame
        Output of utils.multioutput_to_pandas
    col : str, optional
        Column for which to generate a distribution function. The default is
        '[fe/h]'
    range : tuple, optional
        Range in the given column to bin. The default is None, which corresponds
        to the entire range of data
    bin_width : float, optional
        MDF bin width in data units
    """
    if not range:
        range = (stars['[fe/h]'].min(), stars['[fe/h]'].max())
    bins = np.arange(range[0], range[1] + bin_width, bin_width)
    # Calculate remaining stellar mass for each particle
    stars['stellar_mass'] = stars['mass'] * (
        1 - stars['age'].apply(vice.cumulative_return_fraction))
    # Sum remaining stellar mass binned by metallicity
    mdf = stars.groupby([pd.cut(stars['[fe/h]'], bins)])['stellar_mass'].sum()
    mdf /= (mdf.sum() * bin_width)
    return mdf, bins


def box_smooth(hist, bins, width):
    """
    Box-car smoothing function for a pre-generated histogram.

    Parameters
    ----------
    bins : array-like
        Bins dividing the histogram, including the end. Length must be 1 more
        than the length of hist, and bins must be evenly spaced.
    hist : array-like
        Histogram of data
    width : float
        Width of the box-car smoothing function in data units
    """
    bin_width = bins[1] - bins[0]
    box_width = int(width / bin_width)
    box = np.ones(box_width) / box_width
    hist_smooth = np.convolve(hist, box, mode='same')
    return hist_smooth


def smooth_pdf(X, range=None, num=51, kernel='gaussian', bandwidth=0.05):
    """
    Perform a 1D kernel density estimate.

    Parameters
    ----------
    X : array-like
        One-dimensional data
    range : tuple or None, default: None
        Range of data to sample
    num : int, default: 51
        Number of samples to generate
    kernel : str, optional
        Type of kernel to use. The default is 'gaussian'.
    bandwidth : float, optional
        The bandwidth of the kernel. The default is 0.02.

    Returns
    -------
    X_plot : numpy array
        The data coordinates sampled
    dens : numpy array
        The smoothed PDF

    """
    if range == None:
        range = (min(X), max(X))
    X = np.array(X)[:, np.newaxis]
    X_plot = np.linspace(range[0], range[1], num=num,
                         endpoint=True)[:, np.newaxis]
    kde = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(X)
    log_dens = kde.score_samples(X_plot)
    return X_plot[:,0], np.exp(log_dens)


def setup_colorbar(fig, cmap, norm, label=''):
    # fig.subplots_adjust(right=0.95, left=0.12, bottom=0.07, top=0.85,
    #                     wspace=0.05, hspace=0.12)
    # cax = plt.axes([0.12, 0.93, 0.83, 0.02])
    fig.subplots_adjust(right=0.97, left=0.1, bottom=0.2, top=0.93,
                        wspace=0.07, hspace=0.)
    cax = plt.axes([0.15, 0.09, 0.7, 0.02])
    # Add colorbar
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax,
                        orientation='horizontal')
    cbar.set_label(label)#, labelpad=5)
    # cbar.ax.xaxis.set_label_position('top')
    return cax


def discrete_colormap(cmap_name, bounds):
    cmap = plt.get_cmap(cmap_name)
    norm = BoundaryNorm(bounds, cmap.N)
    return cmap, norm


def setup_axes(xlim=FEH_LIM):
    fig, axs = plt.subplots(3, 2, figsize=(3, 4),
                            sharex=True, sharey=True)
    axs[0,0].set_xlim(xlim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.2))
    # Refine axes
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
    return fig, axs


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    migration = sys.argv[1]
    evolution = sys.argv[2]
    RIa = sys.argv[3]
    main(migration, evolution, RIa)
