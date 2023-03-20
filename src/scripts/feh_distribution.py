"""
Plot metallicity distribution functions (MDFs) of [Fe/H] binned by radius.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None # default='warn'
import matplotlib.pyplot as plt
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    import_apogee, apogee_region
from distribution_functions import setup_axes, plot_distributions

MIGRATION = 'diffusion'
FEH_LIM = (-1.1, 0.6)
BIN_WIDTH = 0.01
SMOOTH_WIDTH = 0.2

def main(evolution, RIa, migration=MIGRATION):
    output = '%s/%s/%s' % (migration, evolution, RIa)
    plot_single_comparison(output, verbose=True)


def plot_multiple_comparison(outputs, labels, output_dir=paths.data/'migration',
                             cmap_name='plasma_r', verbose=False, 
                             fname='mdf_feh_multiple.png',
                             double_line_titles=False):
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
        Name of the colormap to use. The default is 'cmap_r'.
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
    fig, axs = setup_axes(ncols=len(outputs)+1, figure_width=7., 
                          xlabel='[Fe/H]', xlim=FEH_LIM, 
                          major_tick_spacing=0.5, cmap_name=cmap_name)
    if double_line_titles:
        # Allow room for two-line axis titles
        fig.subplots_adjust(top=0.9)
    
    # Plot
    for col, output in enumerate(outputs):
        if verbose:
            print('Plotting MDF from %s...' % output)
        stars = multioutput_to_pandas(Path(output_dir) / output)
        plot_distributions(vice_mdf, stars, '[fe/h]', axs[:,col], 
                           label=labels[col], cmap_name=cmap_name)
    
    if verbose:
        print('Plotting MDF from APOGEE...')
    apogee_data = import_apogee()
    plot_distributions(apogee_mdf, apogee_data, 'FE_H', axs[:,col+1], 
                       label='APOGEE DR17', cmap_name=cmap_name)
            
    axs[0,0].set_ylim((0, None))
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
        Name of the colormap to use. The default is 'plasma_r'.
    verbose : bool
        If True, print status updates. The default is False
    fname : str, optional
        The name of the plot output file including its extension. If '', a name
        is automatically generated based on the output directory name. The 
        default is ''.
    """
    # Set up plot
    fig, axs = setup_axes(ncols=2, figure_width=3.25, xlabel='[Fe/H]', 
                          xlim=FEH_LIM, major_tick_spacing=0.5, 
                          cmap_name=cmap_name)
    
    if verbose:
        print('Plotting [Fe/H] distribution from %s...' % output)
    stars = multioutput_to_pandas(Path(output_dir) / output)
    plot_distributions(vice_mdf, stars, '[fe/h]', axs[:,0], 
                       label='VICE', cmap_name=cmap_name)
    
    if verbose:
        print('Plotting APOGEE [Fe/H] distribution...')
    apogee_data = import_apogee()
    plot_distributions(apogee_mdf, apogee_data, 'FE_H', axs[:,1], 
                       label='APOGEE DR17', cmap_name=cmap_name)
    
    # Automatically generate plot filename if none is provided
    if fname == '':
        fname = '%s_%s.png' % tuple(output.split('/')[1:])
            
    axs[0,0].set_ylim((0, None))
    plt.savefig(paths.figures / 'mdf_feh' / fname, dpi=300)
    plt.close()
    if verbose:
        print('Done! Plot is located at src/tex/figures/mdf_feh/%s' % fname)
        

def vice_mdf(stars, col='[fe/h]', galr_lim=(0, 20), absz_lim=(0, 2), 
             xlim=FEH_LIM, bin_width=BIN_WIDTH, smooth_width=SMOOTH_WIDTH):
    """
    Calculate the MDF in [Fe/H] of a region of a VICE multizone output.
    
    Parameters
    ----------
    stars : pandas.DataFrame
        Data from a VICE multizone run.
    col : str, optional
        Column for which to generate a distribution function. The default is
        '[fe/h]'.
    galr_lim : tuple, optional
        Limits of region in galactic radius (kpc). The default is (0, 20).
    absz_lim : tuple, optional
        Limits of region in galactic z-height (kpc). The default is (0, 2).
    xlim : tuple, optional
        Minimum and maximum abundance values. The default is (-1.1, 0.6).
    bin_width : float, optional
        Width of histogram bins in x-axis units. The default is 0.01.
    smooth_width : float, optional
        Width of boxcar smoothing in x-axis units. If 0, the distribution will
        not be smoothed. The default is 0.2.
    
    Returns
    -------
    mdf : numpy.ndarray
        Boxcar-smoothed MDF.
    bin_edges : numpy.ndarray
        [Fe/H] bins including left and right edges, of length len(dndt)+1.
    """
    subset = filter_multioutput_stars(stars, galr_lim, absz_lim, min_mass=0)
    mdf, bin_edges = gen_mdf(subset, col=col, range=xlim, bin_width=bin_width)
    if smooth_width > 0:
        mdf = box_smooth(mdf, bin_edges, smooth_width)
    return mdf, bin_edges


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


def apogee_mdf(data, col='FE_H', galr_lim=(0, 20), absz_lim=(0, 2), 
               xlim=FEH_LIM, bin_width=BIN_WIDTH, smooth_width=SMOOTH_WIDTH):
    """
    Calculate the MDF in [Fe/H] of a region of astroNN data.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Data from APOGEE DR17.
    col : str, optional
        Column name of desired abundance data. The default is 'FE_H'.
    galr_lim : tuple, optional
        Limits of region in galactic radius (kpc). The default is (0, 20).
    absz_lim : tuple, optional
        Limits of region in galactic z-height (kpc). The default is (0, 2).
    xlim : tuple, optional
        Minimum and maximum abundance values. The default is (-1.1, 0.6).
    bin_width : float, optional
        Width of histogram bins in x-axis units. The default is 0.01.
    smooth_width : float, optional
        Width of boxcar smoothing in x-axis units. If 0, the distribution will
        not be smoothed. The default is 0.2.
    
    Returns
    -------
    mdf : numpy.ndarray
        Boxcar-smoothed MDF.
    bin_edges : numpy.ndarray
        [Fe/H] bins including left and right edges, of length len(dndt)+1.
    """
    subset = apogee_region(data, galr_lim, absz_lim)
    bin_edges = np.arange(xlim[0], xlim[1] + bin_width, bin_width)
    mdf, _ = np.histogram(subset[col], bins=bin_edges, density=True)
    if smooth_width > 0:
        mdf = box_smooth(mdf, bin_edges, smooth_width)
    return mdf, bin_edges


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


if __name__ == '__main__':
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main(evolution, RIa)
