"""
This script plots the distribution of stellar ages as predicted by VICE.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars
from apogee_tools import import_apogee, apogee_region
from distribution_functions import setup_axes, plot_distributions
from _globals import DT, END_TIME

MAX_AGE = 14
BIN_WIDTH = 1

def main(evolution, RIa):
    output = 'diffusion/%s/%s' % (evolution, RIa)
    plot_single_comparison(output, verbose=True)
    
    
def plot_multiple_comparison(outputs, labels, output_dir=paths.data/'migration',
                             cmap_name='plasma_r', verbose=False, 
                             fname='adf_multiple_comparison.png',
                             double_line_titles=False):
    """
    One-stop function to generate a complete plot comparing the ADFs of
    multiple VICE multi-zone outputs to the astroNN ADF.
    
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
                          xlabel='Age [Gyr]', xlim=(0, MAX_AGE), 
                          major_tick_spacing=5, cmap_name=cmap_name)
    fig.subplots_adjust(left=0.1)
    if double_line_titles:
        # Allow room for two-line axis titles
        fig.subplots_adjust(top=0.9)
    
    # Plot
    for col, output in enumerate(outputs):
        if verbose:
            print('Plotting age distribution from %s...' % output)
        stars = multioutput_to_pandas(Path(output_dir) / output)
        plot_distributions(vice_adf, stars, axs[:,col], 
                           label=labels[col], cmap_name=cmap_name)
    
    if verbose:
        print('Plotting astroNN age distribution...')
    apogee_data = import_apogee()
    plot_distributions(astroNN_adf, apogee_data, axs[:,-1], 
                       label='astroNN', cmap_name=cmap_name)
            
    axs[0,0].set_ylim((0, None))
    plt.savefig(paths.figures / fname, dpi=300)
    plt.close()
    if verbose:
        print('Done!')
    
    
def plot_single_comparison(output, output_dir=paths.data/'migration',
                           label='VICE', cmap_name='plasma_r', verbose=False,
                           fname=''):
    """
    One-stop function to generate a complete plot comparing the ADF of a single
    VICE output to the astroNN ADF.
    
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
    fig, axs = setup_axes(ncols=2, figure_width=3.25, xlabel='Age [Gyr]', 
                          xlim=(0, MAX_AGE), major_tick_spacing=5, 
                          cmap_name=cmap_name)
    
    if verbose:
        print('Plotting age distribution from %s...' % output)
    stars = multioutput_to_pandas(Path(output_dir) / output)
    plot_distributions(vice_adf, stars, axs[:,0], 
                       label='VICE', cmap_name=cmap_name)
    
    if verbose:
        print('Plotting astroNN age distribution...')
    apogee_data = import_apogee()
    plot_distributions(astroNN_adf, apogee_data, axs[:,1], 
                       label='astroNN', cmap_name=cmap_name)
    
    # Automatically generate plot filename if none is provided
    if fname == '':
        fname = '%s_%s.png' % tuple(output.split('/')[1:])
            
    axs[0,0].set_ylim((0, None))
    plt.savefig(paths.figures / 'adf' / fname, dpi=300)
    plt.close()
    if verbose:
        print('Done! Plot is located at src/tex/figures/%s' % fname)


def vice_adf(stars, galr_lim=(3, 15), absz_lim=(0, 2)):
    """
    Calculate the age distribution function of a region of a VICE multizone
    output.
    
    Parameters
    ----------
    stars : pandas.DataFrame
        Data from a VICE multizone run.
    galr_lim : tuple, optional
        Limits of region in galactic radius (kpc). The default is (3, 15).
    absz_lim : tuple, optional
        Limits of region in galactic z-height (kpc). The default is (0, 2).
    
    Returns
    -------
    dndt : numpy.ndarray
        Distribution of ages.
    bin_edges : numpy.ndarray
        Age bins including left and right edges, of length len(dndt)+1.
    """
    subset = filter_multioutput_stars(stars, galr_lim, absz_lim, min_mass=0)
    dndt, bin_edges = age_distribution(subset, bin_width=BIN_WIDTH)
    return dndt, bin_edges
                
                
def astroNN_adf(data, galr_lim=(3, 15), absz_lim=(0, 2)):
    """
    Calculate the age distribution function of a region of astroNN data.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Data from astroNN (including allStar).
    galr_lim : tuple, optional
        Limits of region in galactic radius (kpc). The default is (3, 15).
    absz_lim : tuple, optional
        Limits of region in galactic z-height (kpc). The default is (0, 2).
    
    Returns
    -------
    dndt : numpy.ndarray
        Distribution of ages.
    bin_edges : numpy.ndarray
        Age bins including left and right edges, of length len(dndt)+1.
    """
    subset = apogee_region(data, galr_lim, absz_lim)
    bin_edges = np.arange(0, MAX_AGE + BIN_WIDTH, BIN_WIDTH)
    dndt, _ = np.histogram(subset['ASTRONN_AGE'], 
                               bins=bin_edges, density=True)
    return dndt, bin_edges


def age_distribution(stars, bin_width=BIN_WIDTH, end_time=END_TIME, dt=DT, 
                     age_col='age', **kwargs):
    """
    Calculate the distribution of ages in a VICE multizone output.
    
    Parameters
    ----------
    stars : DataFrame
        Pandas conversion of a VICE multizone output
    end_time : float, optional
        Simulation end time in Gyr. The default is 13.2
    dt : float, optional
        Simulation time step in Gyr. The default is 0.01
    age_col : str, optional
        Name of column containing ages. The default is 'age'.
    kwargs : dict, optional
        Keyword arguments passed to mean_stellar_mass
        
    Returns
    -------
    dndt : 1xn numpy.ndarray
        Fraction of stars in each age bin
    bins : 1x(n+1) numpy.ndarray
        Age bins derived from simulation timesteps
    """
    bins = np.arange(0, end_time + bin_width, bin_width)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    # Create dummy entries to count at least 0 mass at every age
    temp_df = pd.DataFrame({
        age_col: bins[:-1], 
        'mass': np.zeros(bins[:-1].shape)
    })
    stars = pd.concat([stars, temp_df])
    stars[age_col] = np.round(stars[age_col], decimals=2)
    # Sum stellar mass in each bin
    mass_total, _ = np.histogram(stars[age_col], bins=bins, 
                                 weights=stars['mass'])
    # Calculate remaining stellar mass today
    mass_remaining = mass_total * (1 - np.array(
        [vice.cumulative_return_fraction(age) for age in bin_centers]))
    # Average mass of a star of that particular age
    mass_average = np.array(
        [mean_stellar_mass(age, **kwargs) for age in bin_centers]
    )
    # Number of stars in each age bin
    nstars = np.around(mass_remaining / mass_average)
    # Fraction of stars in each age bin
    if nstars.sum() > 0:
        dndt = nstars / (bin_width * nstars.sum())
    return dndt, bins


def mean_stellar_mass(age, imf=vice.imf.kroupa, mlr=vice.mlr.larson1974,
                      m_lower=0.08, m_upper=100, dm=0.01):
    """
    Calculate the mean mass of a stellar population of a given age.

    Parameters
    ----------
    age : float
        Stellar age in Gyr
    imf : <function>, optional
        Initial mass function which takes mass in solar masses as an argument.
        The default is vice.imf.kroupa
    mlr : <function>, optional
        Mass-lifetime relation which takes age in Gyr as an argument. The
        default is vice.mlr.larson1974
    m_lower : float, optional
        Lower mass limit on IMF in solar masses. The default is 0.08
    m_upper : float, optional
        Upper mass limit on IMF in solar masses. The default is 100
    dm : float, optional
        IMF integration step in solar masses. The default is 0.01

    Returns
    -------
    float
        Mean mass of stars with lifetime greater than or equal to the given age
        weighted by the IMF
    """
    m_max = min((mlr(age, which='age'), m_upper))
    masses = np.arange(m_lower, m_max + dm, dm)
    dndm = np.array([imf(m) for m in masses])
    weighted_mean = np.average(masses, weights=dndm)
    return weighted_mean


if __name__ == '__main__':
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main(evolution, RIa)
