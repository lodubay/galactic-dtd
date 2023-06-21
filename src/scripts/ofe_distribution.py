"""
Plot metallicity distribution functions (MDFs) of [O/Fe] binned by radius.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars
from apogee_tools import import_apogee, apogee_region
from distribution_functions import setup_axes, plot_distributions
from feh_distribution import vice_mdf, apogee_mdf

MIGRATION = 'diffusion'
OFE_LIM = (-0.1, 0.5)
BIN_WIDTH = 0.005
SMOOTH_WIDTH = 0.05


def main(evolution, RIa, migration=MIGRATION):
    output = '%s/%s/%s' % (migration, evolution, RIa)
    plot_single_comparison(output, verbose=True)


def plot_multiple_comparison(outputs, labels, output_dir=paths.data/'migration',
                             cmap_name='plasma_r', verbose=False, 
                             fname='mdf_ofe_multiple.png',
                             double_line_titles=False, figure_width=7.,
                             panel_aspect_ratio=1.5, cbar_width=0.6,
                             smooth_width=SMOOTH_WIDTH):
    """
    One-stop function to generate a complete plot comparing the MDFs of
    multiple VICE multi-zone outputs to the APOGEE MDF.
    
    Parameters
    ----------
    outputs : list
        List of paths to VICE multioutputs starting at the parent directory.
    labels : list
        Titles of each VICE column. Must be same length as outputs.
    output_dir : str or pathlib.Path, optional
        The parent directory for all VICE multizone outputs. The default is
        '../data/migration'.
    cmap_name : str, optional
        Name of the colormap to use. The default is 'plasma_r'.
    verbose : bool, optional
        If True, print status updates. The default is False
    fname : str, optional
        The name of the plot output file including its extension. The default
        is 'adf_multiple_comparison.png'.
    double_line_titles : bool, optional
        If True, provide space for double-line axis titles. The default is 
        False.
    """
        
    # Set up plot
    fig, axs = setup_axes(ncols=len(outputs)+1, figure_width=figure_width, 
                          xlabel='[O/Fe]', xlim=OFE_LIM, 
                          major_tick_spacing=0.2, cmap_name=cmap_name,
                          panel_aspect_ratio=panel_aspect_ratio,
                          cbar_width=cbar_width)
    # Adjust spacing
    fig.subplots_adjust(left=0.08, bottom=0.2)
    if double_line_titles:
        # Allow room for two-line axis titles
        fig.subplots_adjust(top=0.9)
    
    # Plot
    if verbose:
        print('Plotting [O/Fe] distribution from APOGEE...')
    apogee_data = import_apogee()
    plot_distributions(apogee_mdf, apogee_data, 'O_FE', axs[:,0], 
                       label='APOGEE DR17', cmap_name=cmap_name,
                       func_kwargs={'xlim': OFE_LIM, 'bin_width': BIN_WIDTH,
                                    'smooth_width': smooth_width})
    
    for col, output in enumerate(outputs):
        if verbose:
            print('Plotting [O/Fe] distribution from %s...' % output)
        stars = multioutput_to_pandas(Path(output_dir) / output)
        plot_distributions(vice_mdf, stars, '[o/fe]', axs[:,col+1], 
                           label=labels[col], cmap_name=cmap_name,
                           func_kwargs={'xlim': OFE_LIM, 'bin_width': BIN_WIDTH,
                                        'smooth_width': smooth_width})
            
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    plt.savefig(paths.figures / fname, dpi=300)
    plt.close()
    if verbose:
        print('Done!')
    
    
def plot_single_comparison(output, output_dir=paths.data/'migration',
                           label='VICE', cmap_name='plasma_r', verbose=False,
                           fname=''):
    """
    One-stop function to generate a complete plot comparing the MDF of a single
    VICE output to the APOGEE MDF.
    
    Parameters
    ----------
    output : str
        Path to VICE multioutput starting at the parent directory.
    output_dir : str or pathlib.Path, optional
        The parent directory for all VICE multizone outputs. The default is
        '../data/migration'.
    label : str, optional
        Title of VICE column. The default is 'VICE'.
    cmap_name : str
        Name of the colormap to use. The default is 'cmap_r'.
    verbose : bool
        If True, print status updates. The default is False
    fname : str, optional
        The name of the plot output file including its extension. If '', a name
        is automatically generated based on the output directory name. The 
        default is ''.
    """
    # Set up plot
    fig, axs = setup_axes(ncols=2, figure_width=3.25, xlabel='[O/Fe]', 
                          xlim=OFE_LIM, major_tick_spacing=0.2, 
                          cmap_name=cmap_name)
    
    if verbose:
        print('Plotting [O/Fe] distribution from %s...' % output)
    stars = multioutput_to_pandas(Path(output_dir) / output)
    plot_distributions(vice_mdf, stars, '[o/fe]', axs[:,0], 
                       label='VICE', cmap_name=cmap_name,
                       func_kwargs={'xlim': OFE_LIM, 'bin_width': BIN_WIDTH,
                                    'smooth_width': SMOOTH_WIDTH})
    
    if verbose:
        print('Plotting APOGEE [O/Fe] distribution...')
    apogee_data = import_apogee()
    plot_distributions(apogee_mdf, apogee_data, 'O_FE', axs[:,1], 
                       label='APOGEE DR17', cmap_name=cmap_name,
                       func_kwargs={'xlim': OFE_LIM, 'bin_width': BIN_WIDTH,
                                    'smooth_width': SMOOTH_WIDTH})
    
    # Automatically generate plot filename if none is provided
    if fname == '':
        fname = '%s_%s.png' % tuple(output.split('/')[1:])
            
    for ax in axs[:,0]:
        ax.set_ylim((0, None))
    figpath = paths.debug / 'ofe_df' / output.split('/')[0] / fname
    plt.savefig(figpath, dpi=300)
    plt.close()
    if verbose:
        print('Done! Plot is located at ', str(figpath))


if __name__ == '__main__':
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main(evolution, RIa)
