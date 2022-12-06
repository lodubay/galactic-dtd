"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights with
overlayed contours of APOGEE abundances. Takes a couple minutes to run.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import multioutput_to_pandas
from age_ofe_vice import plot_age_ofe_stars, plot_post_process_tracks, \
    setup_axes
from age_ofe_apogee import plot_medians, plot_contours
from utils import import_astroNN, select_giants
import paths

def main(output_name, migration_dir='../data/migration',
         stars_cmap='winter', apogee_cmap='magma'):
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
    print('Importing VICE')
    stars = multioutput_to_pandas(output_name, migration_dir)
    # Plot simulation output
    print('Plotting VICE stars')
    fig, axs = plot_age_ofe_stars(stars, stars_cmap)
    # Import APOGEE data
    print('Importing APOGEE')
    apogee_data = select_giants(import_astroNN())
    apogee_data.dropna(subset=['O_FE', 'ASTRONN_AGE'], inplace=True)
    # print('Plotting APOGEE contours')
    # plot_contours(axs, apogee_data, apogee_cmap, linewidths=0.5)
    print('Plotting APOGEE medians')
    plot_medians(axs, apogee_data)
    # Add post-process abundance track
    print('Plotting abundance tracks')
    # plot_post_process_track(output_name, axs, galr=8, data_dir=migration_dir)
    plot_post_process_tracks(output_name, axs, data_dir=migration_dir)
    # fig.suptitle(output_name)
    migration, evolution, RIa = output_name.split('/')[-3:]
    plt.savefig(paths.figures / ('age_ofe/%s_%s.png' % (evolution, RIa)), dpi=300)
    # plt.savefig(paths.figures / ('ofe_feh_%s_%s.pdf' % (evolution, RIa)), dpi=300)
    plt.close()


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))
