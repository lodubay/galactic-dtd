"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights with
additional contours of APOGEE abundances.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import multioutput_to_pandas
from ofe_feh_vice import plot_ofe_feh_stars, plot_post_process_track
from ofe_feh_vice import GALR_BINS, ABSZ_BINS
from ofe_feh_apogee import kde2D, apogee_region

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
    # Import APOGEE data
    apogee_data = pd.read_csv(Path(apogee_path))
    plot_apogee_contours(axs, apogee_data, apogee_cmap)
    # Add post-process abundance track
    plot_post_process_track(output_name, axs, galr=8, data_dir=migration_dir)
    fig.suptitle(output_name)
    plt.savefig('ofe_feh_compare_%s.pdf' % output_name.split('/')[-1], dpi=300)
    plt.close()


def plot_apogee_contours(axs, data, cmap='magma'):
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
            subset = apogee_region(data, galr_lim, absz_lim)
            x, y, logz = kde2D(subset['FE_H'], subset['O_FE'], 0.03)
            # Scale by total number of stars in region
            logz += np.log(subset.shape[0])
            # Contour levels in log-likelihood space
            levels = np.arange(10, 14, 0.5)
            ax.contour(x, y, logz, levels, cmap=cmap)


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))