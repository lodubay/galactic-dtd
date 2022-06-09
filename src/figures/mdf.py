"""
Plot metallicity distribution functions (MDFs) of [O/Fe] binned by radius.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import vice
from ofe_feh_vice import GALR_BINS, ABSZ_BINS, FEH_LIM, OFE_LIM, ZONE_WIDTH
from ofe_feh_apogee import apogee_region
from utils import multioutput_to_pandas, filter_multioutput_stars

def main(output_name, migration_dir='../data/migration_outputs',
         apogee_path='../data/APOGEE/dr17_cut_data.csv', cmap_name='copper'):
    print('Importing data')
    stars = multioutput_to_pandas(output_name, migration_dir)
    data = pd.read_csv(Path(apogee_path))

    xlim = (-0.7, 0.7)
    x_plot = np.linspace(xlim[0], xlim[1], 101, endpoint=True)[:, np.newaxis]
    ylim = (0, 4.5)
    fig, axs = plt.subplots(3, 2, figsize=(6, 9), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    cmap = plt.get_cmap(cmap_name)
    colors = cmap([r/15 for r in GALR_BINS[:-1]])
    for i in range(len(ABSZ_BINS)-1):
        print('Panel %s' % (i+1))
        absz_lim = ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i]
        axs[i,0].text(0.55, 0.85, r'$|z| = %s - %s$' % tuple(absz_lim),
                      transform=axs[i,0].transAxes)
        for j in range(len(GALR_BINS)-1):
            galr_lim = GALR_BINS[j:j+2]
            label = '%s - %s kpc' % tuple(galr_lim)
            # Plot VICE in left panels
            vice_subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
                                                   ZONE_WIDTH, min_mass=0)
            X_plot, vice_pdf = smooth_pdf(vice_subset['[fe/h]'], range=xlim)
            axs[i,0].plot(X_plot, vice_pdf, color=colors[j])
            # axs[i,0].hist(vice_subset['[o/h]'], bins=80, range=xlim,
            #             density=True, histtype='step', color=colors[j])
            # mdf_sum = 0
            # for z in range(int(galr_lim[0] / ZONE_WIDTH), int(galr_lim[1] / ZONE_WIDTH)):
            #     singlezone_output = Path('%s.vice/zone%s' % (output_name, z))
            #     mdf = vice.mdf(str(Path(migration_dir) / singlezone_output))
            #     mdf_sum += np.array(mdf['dn/d[o/fe]'])
            # bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
            # axs[i,0].hist(bins[:-1], bins, weights=mdf_sum, histtype='step',
            #               color=colors[j], density=True)

            # Plot APOGEE in right panels
            apogee_subset = apogee_region(data, galr_lim, absz_lim)
            X_plot, apogee_pdf = smooth_pdf(apogee_subset['FE_H'], range=xlim)
            # x = np.array(apogee_subset['O_FE'] + apogee_subset['FE_H'])
            # x = x[:, np.newaxis]
            # kde = KernelDensity(kernel='gaussian', bandwidth=0.02).fit(x)
            # apogee_pdf = np.exp(kde.score_samples(x_plot))
            axs[i,1].plot(X_plot, apogee_pdf, color=colors[j])
            # axs[i,1].hist(apogee_subset['O_FE'] + apogee_subset['FE_H'], bins=80, range=xlim,
            #             density=True, histtype='step', color=colors[j])
    axs[0,0].set_xlim(xlim)
    axs[0,0].set_ylim(ylim)
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    for ax in axs[:,0]:
        ax.set_ylabel('PDF')
    axs[0,0].set_title(output_name)
    axs[0,1].set_title('APOGEE DR17')
    plt.savefig('feh_mdf.png', dpi=300)
    plt.close()


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

if __name__ == '__main__':
    main('diffusion/insideout/long_delay')