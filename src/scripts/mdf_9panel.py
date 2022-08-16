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
from sklearn.neighbors import KernelDensity
import vice
from _globals import ABSZ_BINS, ZONE_WIDTH, GALR_BINS
from ofe_feh_apogee import apogee_region
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    import_allStar, discrete_colormap, setup_discrete_colorbar
import paths

global FEH_LIM
global BIN_WIDTH
global SMOOTH_WIDTH

FEH_LIM = (-1.1, 0.6)
BIN_WIDTH = 0.01
SMOOTH_WIDTH = 0.2

def main(evolution, RIa, cmap_name='plasma_r'):
    output = '%s/%s' % (evolution, RIa)
    stars_pp = multioutput_to_pandas(output, paths.data / 'migration' / 'post-process')
    stars_diff = multioutput_to_pandas(output, paths.data / 'migration' / 'diffusion')
    apogee_data = import_allStar()

    fig, axs = setup_axes(xlim=FEH_LIM)
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_discrete_colorbar(fig, cmap, norm, 
                                  label=r'Galactocentric radius [kpc]')
    # Define color scheme
    rmin, rmax = GALR_BINS[0], GALR_BINS[-2]
    colors = cmap([(r-rmin)/(rmax-rmin) for r in GALR_BINS[:-1]])

    for i in range(len(ABSZ_BINS)-1):
        absz_lim = ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i]
        axs[i,0].set_ylabel(r'$|z| = %s - %s$' % tuple(absz_lim))
        for j in range(len(GALR_BINS)-1):
            galr_lim = GALR_BINS[j:j+2]
            # Plot VICE post-process in left panels
            pp_subset = filter_multioutput_stars(stars_pp,
                                                 galr_lim, absz_lim,
                                                 ZONE_WIDTH, min_mass=0)
            pp_mdf, bins = gen_mdf(pp_subset, range=FEH_LIM,
                                   bin_width=BIN_WIDTH)
            pp_smooth = box_smooth(pp_mdf, bins, SMOOTH_WIDTH)
            x_plot = bins[:-1] + BIN_WIDTH/2
            axs[i,0].plot(x_plot, pp_smooth, color=colors[j], linewidth=1)

            # Plot VICE diffusion in center panels
            diff_subset = filter_multioutput_stars(stars_diff,
                                                   galr_lim, absz_lim,
                                                   ZONE_WIDTH, min_mass=0)
            diff_mdf, bins = gen_mdf(diff_subset, range=FEH_LIM,
                                     bin_width=BIN_WIDTH)
            diff_smooth = box_smooth(diff_mdf, bins, SMOOTH_WIDTH)
            axs[i,1].plot(x_plot, diff_smooth, color=colors[j], linewidth=1)

            # Plot APOGEE in right panels
            apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
            apogee_mdf, _ = np.histogram(apogee_subset['FE_H'], bins=bins,
                                         density=True)
            apogee_smooth = box_smooth(apogee_mdf, bins, SMOOTH_WIDTH)
            axs[i,2].plot(x_plot, apogee_smooth, color=colors[j], linewidth=1)

    ylim = axs[0,0].get_ylim()
    axs[0,0].set_ylim((1.5*ylim[0], None))
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    # for ax in axs[:,0]:
        # ax.set_ylabel('PDF')
    # axs[0,0].set_title('/'.join((evolution, RIa)))
    # axs[0,0].set_title('VICE')
    axs[0,0].set_title('Post-process')
    axs[0,1].set_title('Diffusion')
    axs[0,2].set_title('APOGEE DR17')
    plt.savefig(paths.figures / ('mdf_%s_%s.png' % (evolution, RIa)), dpi=300)
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
        range = (stars[col].min(), stars[col].max())
    bins = np.arange(range[0], range[1] + bin_width, bin_width)
    # Calculate remaining stellar mass for each particle
    stars['stellar_mass'] = stars['mass'] * (
        1 - stars['age'].apply(vice.cumulative_return_fraction))
    # Sum remaining stellar mass binned by metallicity
    mdf = stars.groupby([pd.cut(stars[col], bins)])['stellar_mass'].sum()
    # Normalize
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


def setup_axes(xlim=FEH_LIM, ncols=3, absz_bins=ABSZ_BINS, galr_bins=GALR_BINS):
    nrows = len(absz_bins)
    fig, axs = plt.subplots(3, ncols, figsize=(1.5*ncols, nrows+1),
                            sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, top=0.93, right=0.97, bottom=0.1,
                        wspace=0.07, hspace=0.)
    # axs[0,0].set_xlim((xlim[0]-0.09, xlim[1]+0.09))
    axs[0,0].set_xlim(xlim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # axs[0,0].yaxis.set_major_locator(MultipleLocator(1))
    # axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.2))
    # Refine axes
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_bounds(xlim[0], xlim[1])
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
    # Turn on y-axis spine and ticks for left-hand panels
    # for ax in axs[:,0]:
    #     ax.spines['left'].set_visible(True)
    #     ax.yaxis.set_ticks_position('left')
    #     ax.spines['left'].set_bounds(0, 2.6)
    return fig, axs


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main(evolution, RIa)
