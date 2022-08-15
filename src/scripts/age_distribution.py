"""
This script plots the distribution of stellar ages as predicted by VICE.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars, import_astroNN, get_bin_centers
from _globals import DT, GALR_BINS, ZONE_WIDTH, ABSZ_BINS, END_TIME
from mdf_9panel import setup_colorbar, discrete_colormap, box_smooth
from ofe_feh_apogee import apogee_region

MAX_AGE = 14
BIN_WIDTH = 1
# SMOOTH_WIDTH = 1

def main(evolution, RIa):
    output = 'diffusion/%s/%s' % (evolution, RIa)
    plot_single_comparison(output)
    
    
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
    cmap_name : str
        Name of the colormap to use. The default is 'cmap_r'.
    verbose : bool
        If True, print status updates. The default is False
    fname : str, optional
        The name of the plot output file including its extension. The default
        is 'adf_multiple_comparison.png'.
    double_line_titles : bool, optional
        If True, provide space for double-line axis titles. The default is 
        False.
    """
        
    # Set up plot
    fig, axs = setup_axes(ncols=len(outputs)+1, figure_width=7.)
    if double_line_titles:
        # Allow room for two-line axis titles
        fig.subplots_adjust(top=0.9)
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_colorbar(fig, cmap, norm, label=r'Galactocentric radius [kpc]')
    # Define color scheme
    colors = get_color_list(cmap, GALR_BINS)
    # Plot
    for col, output in enumerate(outputs):
        if verbose:
            print('Plotting age distribution from %s...' % output)
        stars = multioutput_to_pandas(Path(output_dir) / output)
        plot_vice_adf(stars, axs[:,col], colors=colors, label=labels[col])
    
    if verbose:
        print('Plotting astroNN age distribution...')
    astroNN_data = import_astroNN()
    plot_astroNN_adf(astroNN_data, axs[:,col+1], colors=colors)
            
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
    fig, axs = setup_axes(ncols=2, figure_width=7.)
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_colorbar(fig, cmap, norm, label=r'Galactocentric radius [kpc]')
    colors = get_color_list(cmap, GALR_BINS)
    
    if verbose:
        print('Plotting age distribution from %s...' % output)
    stars = multioutput_to_pandas(Path(output_dir) / output)
    plot_vice_adf(stars, axs[:,0], colors=colors, label=label)
    
    if verbose:
        print('Plotting astroNN age distribution...')
    astroNN_data = import_astroNN()
    plot_astroNN_adf(astroNN_data, axs[:,1], colors=colors)
    
    # Automatically generate plot filename if none is provided
    if fname == '':
        fname = 'adf_%s_%s.png' % tuple(output.split('/')[1:])
            
    axs[0,0].set_ylim((0, None))
    plt.savefig(paths.figures / fname, dpi=300)
    plt.close()
    if verbose:
        print('Done! Plot is located at src/tex/figures/%s' % fname)    
    

def plot_astroNN_adf(data, axs, colors=[], label='astroNN', 
                     age_bin_width=BIN_WIDTH, max_age=MAX_AGE,
                     absz_bins=ABSZ_BINS, galr_bins=GALR_BINS):
    """
    Plot age distributions from astroNN data.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Combined data from APOGEE and astroNN.
    axs : list of matplotlib.axes.Axes
        Axes on which to plot the age distributions, the length of which must
        correspond to len(absz_bins)-1 (3 by default); usually a single
        column from a larger array of axes.
    colors : list, optional
        List of colors corresponding to galactic radius. If len(colors) !=
        len(galr_bins)-1, the default matplotlib color scheme will be used.
        The default is [].
    label : str, optional
        Axis column label. The default is 'astroNN'.
    age_bin_width : float, optional
        Width of age bins in Gyr. The default is 1.
    max_age : float, optional
        Maximum age in Gyr to include. The default is 14.
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
    galr_bins : list, optional
        Bin edges of galactic radius in kpc. The default is
        [3, 5, 7, 9, 11, 13, 15].
    """
    if len(colors) != len(galr_bins) - 1:
        colors = [None for galr in galr_bins[:-1]]
    if len(axs) == len(absz_bins) - 1:
        for i, ax in enumerate(axs.flatten()):
            absz_lim = absz_bins[-(i+2):len(absz_bins)-i]
            for j in range(len(galr_bins)-1):
                galr_lim = galr_bins[j:j+2]
                subset = apogee_region(data, galr_lim, absz_lim)
                age_bins = np.arange(0, max_age + age_bin_width, age_bin_width)
                age_hist, _ = np.histogram(subset['ASTRONN_AGE'], 
                                           bins=age_bins, density=True)
                ax.plot(get_bin_centers(age_bins), age_hist, 
                        color=colors[j], linewidth=1)
        axs[0].set_title(label)
    else:
        raise ValueError('Mismatch between axes and z-height bins.')
    
    
def plot_vice_adf(stars, axs, colors=[], label='VICE', age_bin_width=BIN_WIDTH,
                  absz_bins=ABSZ_BINS, galr_bins=GALR_BINS, 
                  zone_width=ZONE_WIDTH):
    """
    Plot age distributions from a VICE multi-zone run.
    
    Parameters
    ----------
    stars : pandas.DataFrame
        Output of utils.multioutput_to_stars.
    axs : list of matplotlib.axes.Axes
        Axes on which to plot the age distributions, the length of which must
        correspond to len(absz_bins)-1 (3 by default); usually a single
        column from a larger array of axes.
    colors : list, optional
        List of colors corresponding to galactic radius. If len(colors) !=
        len(galr_bins)-1, the default matplotlib color scheme will be used.
        The default is [].
    label : str, optional
        Axis column label. The default is 'VICE'.
    age_bin_width : float, optional
        Width of age bins in Gyr. The default is 1.
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
    galr_bins : list, optional
        Bin edges of galactic radius in kpc. The default is 
        [3, 5, 7, 9, 11, 13, 15].
    zone_width : float, optional
        Width of VICE zones in kpc. The default is 0.1.
    """
    if len(colors) != len(galr_bins) - 1:
        colors = [None for galr in galr_bins[:-1]]
    if len(axs) == len(absz_bins) - 1:
        for i, ax in enumerate(axs.flatten()):
            absz_lim = absz_bins[-(i+2):len(absz_bins)-i]
            for j in range(len(galr_bins)-1):
                galr_lim = galr_bins[j:j+2]
                subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
                                                  zone_width, min_mass=0)
                dndt, bins = age_distribution(subset, bin_width=age_bin_width)
                ax.plot(get_bin_centers(bins), dndt, 
                        color=colors[j], linewidth=1)
        axs[0].set_title(label)
    else:
        raise ValueError('Mismatch between axes and z-height bins.')

def get_color_list(cmap, bins):
    """
    Split a discrete colormap into a list of colors based on bin edges.
    
    Parameters
    ----------
    cmap : matplotlib colormap
    bins : array-like
        Bin edges, including left- and right-most edges
    
    Returns
    -------
    list
        List of colors of length len(bins) - 1
    """
    rmin, rmax = bins[0], bins[-2]
    colors = cmap([(r-rmin)/(rmax-rmin) for r in bins[:-1]])
    return colors

def setup_axes(ncols=2, figure_width=3.25, xlim=(0, MAX_AGE), 
               absz_bins=ABSZ_BINS):
    """
    Set up matplotlib figure and axes for the age distribution plot.
    
    Parameters
    ----------
    ncols : int, optional
        Number of columns in figure. The default is 2.
    figure_width : float, optional
        Width of the figure in inches. For an AASTeX document, this should be
        3.25 for a single-column figure or 7.0 for a double-column figure. The
        figure height will be scaled to maintain a consistent aspect ratio.
        The default is 3.25.
    xlim : tuple or list, optional
        Limits of x-axis in Gyr. The default is (0, 14).
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : list of matplotlib.axes.Axes
    """
    # Determine figure dimensions
    nrows = len(absz_bins) - 1
    ax_width = (figure_width - 0.25) / ncols
    ax_height = ax_width / 1.5
    figure_height = ax_height * nrows + 1
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True,
                            figsize=(figure_width, figure_height))
    fig.subplots_adjust(left=0.07, top=0.93, right=0.97, bottom=0.1,
                        wspace=0.07, hspace=0.)
    # Format x-axis
    axs[0,0].set_xlim(xlim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    for ax in axs[-1]:
        ax.set_xlabel('Age [Gyr]')
    # Remove spines and y-axis labels
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
    # Label rows
    for i in range(len(absz_bins)-1):
        absz_lim = tuple(absz_bins[-(i+2):len(absz_bins)-i])
        axs[i,0].set_ylabel(r'$|z| = %s - %s$' % absz_lim)
    return fig, axs


def age_distribution(stars, bin_width=BIN_WIDTH, end_time=END_TIME, dt=DT, 
                     **kwargs):
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
    temp_df = pd.DataFrame({'age': bins[:-1], 'mass': np.zeros(bins[:-1].shape)})
    stars = pd.concat([stars, temp_df])
    stars['age'] = np.round(stars['age'], decimals=2)
    # Sum stellar mass in each bin
    mass_total, _ = np.histogram(stars['age'], bins=bins, weights=stars['mass'])
    # mass_total = stars.groupby(['age']).sum()['mass']
    # ages = np.array(mass_total.index)
    # mass_total = np.array(mass_total)
    # Calculate remaining stellar mass today
    mass_remaining = mass_total * (1 - np.array(
        [vice.cumulative_return_fraction(age) for age in bin_centers]))
    # Average mass of a star of that particular age
    mass_average = np.array([mean_stellar_mass(age, **kwargs) for age in bin_centers])
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


# def f_survive(age, mlr='larson1974', imf='kroupa', m_lower=0.08, m_upper=100,
#               dm=0.01):
#     """
#     Calculate the surviving mass fraction of a stellar population with a given
#     age.

#     Parameters
#     ----------
#     age : float
#         Age of the stellar population in Gyr
#     mlr : str
#         Mass-lifetime relation (MLR) to use. The default is 'larson1974'.
#     imf : str
#         Which IMF to use. Options are 'kroupa' or 'salpeter'.
#     m_lower : float
#         Lower limit of stellar mass. The default is 0.08 solar masses.
#     m_upper : float
#         Upper limit of stellar mass. The default is 100 solar masses.
#     dm : float
#         Integration step size in solar masses. The default is 0.01.

#     Returns
#     -------
#     float
#         Stellar population surviving mass fraction.
#     """
#     mlr_select = {
#         'larson1974': vice.mlr.larson1974,
#         'mm1989': vice.mm1989,
#         'pm1993': vice.mlr.pm1993,
#         'ka1997': vice.mlr.ka1997,
#         'hpt2000': vice.mlr.hpt2000,
#         'vincenzo2016': vice.mlr.vincenzo2016,
#         'powerlaw': vice.mlr.powerlaw
#     }
#     if mlr in mlr_select.keys():
#         # Mass of a star with a lifetime equal to age
#         mass = mlr_select[mlr](age, which='age')
#         m_arr = np.array(m_lower, min((m_upper, mass)), dm)
#         normal_imf = NormalIMF(which=imf, m_lower=m_lower, m_upper=m_upper,
#                                dm=dm)
#         f_survive = sum([normal_imf(m) for m in m_arr])
#         return f_survive
#     else:
#         raise ValueError('MLR must be in acceptable list.')


if __name__ == '__main__':
    main('insideout_johnson21', 'powerlaw_slope11_delay040')
