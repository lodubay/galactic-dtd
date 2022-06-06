"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights with
additional contours of APOGEE abundances.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from sklearn.neighbors import KernelDensity
from utils import multioutput_to_pandas
from ofe_feh_vice import plot_ofe_feh_stars, plot_post_process_track
from ofe_feh_vice import GALR_BINS, ABSZ_BINS

def main(output_name, migration_dir='../data/migration_outputs',
         stars_cmap='winter', apogee_path='../data/APOGEE/dr17_cut_data.csv',
         apogee_cmap='magma'):
    """
    Parameters
    ----------
    output_name : str
        Path to the VICE multizone output
    migration_dir : str, optional
        Directory containing all VICE multizone outputs
    stars_cmap : str, optional [default: 'winter']
        Name of colormap to use for VICE stars output
    apogee_path : str, optional
        Path to APOGEE data
    apogee_cmap : str, optional [default: 'magma']
        Name of colormap to use for APOGEE data contours
    """
    # Import multioutput stars data
    stars = multioutput_to_pandas(output_name, migration_dir)
    # Plot simulation output
    fig, axs = plot_ofe_feh_stars(stars, stars_cmap)
    axs = plot_apogee_contours(axs, apogee_path, apogee_cmap)
    plot_post_process_track(output_name, axs, galr=8, data_dir=migration_dir)
    plt.savefig('ofe_feh_compare.pdf', dpi=300)
    plt.close()


def plot_apogee_contours(axs, filename, cmap_name='magma'):
    """
    Add contours of APOGEE abundances on top of VICE scatterplot.

    Parameters
    ----------
    axs : list of axes
    filename : str
        Path to APOGEE data file
    cmap_name : str
        Name of colormap to apply to contours
    """
    cmap = get_cmap(cmap_name)
    data = pd.read_csv(Path(filename))
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = apogee_region(data, galr_lim, absz_lim)
            x, y, z = kde2D(subset['FE_H'], subset['O_FE'], 0.05)
            ax.contour(x, y, z, cmap=cmap)
    return axs


def apogee_region(data, galr_lim=(0, 20), absz_lim=(0, 5)):
    """
    Slice astroNN data within a given Galactic region of radius and z-height.

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


def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs):
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


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))