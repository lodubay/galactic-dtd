"""
This script plots a grid of panels showing [O/Fe] vs age in multiple galactic
regions. The plot includes a sample of VICE stellar populations, mass-weighted 
medians of the VICE output, and the median estimated ages from Feuillet et al. 
(2019), Mackereth et al. (2019; astroNN), or Leung et al. (2023).
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# import vice
from scatter_plot_grid import setup_axes, setup_colorbar
from utils import weighted_quantile, group_by_bins, feuillet2019_data
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee, apogee_region
from _globals import GALR_BINS, ABSZ_BINS
import paths

GALR_BINS = GALR_BINS[:-1]
AGE_LIM_LINEAR = (-1, 15)
AGE_LIM_LOG = (0.3, 20)
OFE_LIM = (-0.15, 0.55)
OFE_BIN_WIDTH = 0.05
AGE_SOURCES = ['F19', # Feuillet et al. 2019, solar neighborhood only
               'M19', # Mackereth et al. 2019, astroNN
               'L23'] # Leung et al. 2023, variational encoder-decoder
AGE_LABELS = {'F19': 'Feuillet et al. (2019)',
              'M19': 'Mackereth et al. (2019)',
              'L23': 'Leung et al. (2023)'}

def main(evolution, RIa, migration='gaussian', verbose=False, **kwargs):
    # Import VICE multi-zone output data
    output_name = '/'.join([migration, evolution, RIa, 'diskmodel'])
    # Import APOGEE and astroNN data
    apogee_data = import_apogee(verbose=verbose)
    # Main plot function
    plot_age_ofe(output_name, apogee_data, verbose=verbose, **kwargs)


def plot_age_ofe(output_name, apogee_data, fname='age_ofe.png', ages='L23', 
                 cmap='winter', log=True, score=False, uncertainties=True, 
                 verbose=False):
    """
    Plot a grid of [O/Fe] vs age across multiple Galactic regions.
    
    Parameters
    ----------
    output_name : str
        Name of multizone output within the data_dir directory.
    apogee_data : pandas.DataFrame
    fname : str, optional
        File name (excluding parent directory) of plot output. The default is
        'age_ofe.png'.
    ages : str, optional
        Source for age data. Must be one of ('F19', 'M19', 'L23'). The default
        is 'L23'.
    cmap : str, optional
        Name of colormap for scatterplot. The default is 'winter'.
    log : bool, optional
        If True, plot the x (age) axis on a log scale. The default is False.
    score : bool, optional
        If True, calculate a numerical score for how well the VICE output
        matches the age data. The default is False.
    uncertainties : bool, optional
        If True, convolve VICE output data with observational uncertainties
        from APOGEE and age data. The default is False.
    verbose : bool, optional
        If True, print verbose output to the terminal. The default is False.
    savedir : str or pathlib.Path, optional
        Parent directory to save the plot. The default is ../debug/age_ofe/.
        
    Returns
    -------
    float
        If scores==True, the weighted average score across all regions. The
        score represents the difference between the median ages from VICE
        and APOGEE. A lower score represents a better fit.
    """
    # Error handling
    if ages not in AGE_SOURCES:
        raise ValueError('Parameter "ages" must be in %s.' % AGE_SOURCES)
    # Import VICE multizone outputs
    vice_stars = MultizoneStars.from_output(output_name, verbose=verbose)
    # Model uncertainties
    if uncertainties:
        vice_stars.model_uncertainty(apogee_data, age_source=ages, inplace=True)
    # Age data source
    if ages == 'L23':
        age_col = 'LATENT_AGE'
    else:
        age_col = 'ASTRONN_AGE'
    
    # Feuillet+ 2019 ages require specific regions
    if ages == 'F19':
        local_galr_bins = [5, 7, 9, 11, 13]
        local_absz_bins = [0, 0.5, 1, 2]
    else:
        local_galr_bins = GALR_BINS
        local_absz_bins = ABSZ_BINS
    
    # Set x-axis limits
    if log:
        age_lim = AGE_LIM_LOG
    else:
        age_lim = AGE_LIM_LINEAR
    
    # Set up figure
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig, axs = setup_axes(xlim=age_lim, ylim=OFE_LIM, 
                          xlabel='Age [Gyr]', ylabel='[O/Fe]',
                          xlabelpad=2, ylabelpad=4,
                          galr_bins=local_galr_bins, absz_bins=local_absz_bins)
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]', pad=0.02)
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    scores = []
    weights = []
    
    # Plot sampled points and medians
    if verbose: 
        print('Plotting [O/Fe] vs age in galactic regions...')
    for i, row in enumerate(axs):
        absz_lim = (local_absz_bins[-(i+2)], local_absz_bins[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (local_galr_bins[j], local_galr_bins[j+1])
            if verbose:
                print('\tRGal=%s kpc, |z|=%s kpc' \
                      % (str(galr_lim), str(absz_lim)))
            vice_subset = vice_stars.region(galr_lim, absz_lim)
            vice_subset.scatter_plot(ax, 'age', '[o/fe]', color='galr_origin',
                                     cmap=cmap, norm=cbar.norm)
            # Plot Feuillet+ 2019 ages
            if ages == 'F19':
                plot_feuillet2019(ax, galr_lim, absz_lim)
            else:
                apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
                plot_astroNN_medians(ax, apogee_subset, 
                                     age_col=age_col, label=AGE_LABELS[ages])
            plot_vice_medians(ax, vice_subset.stars, label='Model')
            # Score how well the distributions align based on the RMS of the
            # difference in medians in each [O/Fe] bin
            if score and (ages in ('M19', 'L23')):
                d_rms = rms_median_diff(vice_subset.stars, apogee_subset, 
                                        age_col=age_col)
                ax.text(0.07, 0.67, r'$\Delta_{\rm{RMS}}=%.02f$' % d_rms, 
                        transform=ax.transAxes)
                scores.append(d_rms)
                weights.append(apogee_subset.shape[0])
            # Add legend to top-right panel
            if i==0 and j==len(row)-1:
                ax.legend(loc='upper left', frameon=False, handletextpad=0.1,
                          borderpad=0.2, handlelength=1.5)
                             
    # Set x-axis scale and ticks
    if log:
        axs[0,0].set_xscale('log')
        axs[0,0].xaxis.set_major_formatter(FormatStrFormatter('%d'))
    else:
        axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
        
    # Set y-axis ticks
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    # Output figure
    if uncertainties:
        fname = 'age_ofe_%s_errors.png' % ages
    else:
        fname = 'age_ofe_%s.png' % ages
    save_dir = paths.figures / 'supplementary' / output_name.split('/diskmodel')[0]
    if not save_dir.exists():
        save_dir.mkdir(parents=True)
    fig.savefig(save_dir / fname, dpi=300)
    plt.close()
    
    # Return overall model score
    if score:
        return np.average(scores, weights=weights)


def rms_median_diff(vice_stars, apogee_data, age_col='LATENT_AGE', 
                    ofe_lim=OFE_LIM, ofe_bin_width=OFE_BIN_WIDTH):
    """
    Calculate the RMS of the difference in medians between VICE and data ages.
    
    Parameters
    ----------
    vice_stars : pandas.DataFrame
        VICE multizone stars data, typically a subset of a galactic region.
    apogee_data : pandas.DataFrame
        APOGEE data with ages, typically a subset of a galactic region.
    age_col : str, optional
        Name of column with age data in apogee_data. The default is 'LATENT_AGE'
        which is the Leung et al. (2023) ages.
    ofe_lim : tuple, optional
        Outermost bounds on [O/Fe]. The default is (-0.15, 0.55).
    ofe_bin_width : float, optional
        The [O/Fe] bin width in dex. The default is 0.05.
        
    Returns
    -------
    d_rms : float
        Root-mean-square of the difference in median ages in each [O/Fe] bin.
    """
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    # bin APOGEE ages by [O/Fe]
    apogee_grouped = group_by_bins(apogee_data, 'O_FE', ofe_bins)[age_col]
    apogee_medians = apogee_grouped.median()
    # count all APOGEE stars in each bin
    apogee_counts = apogee_grouped.count()
    # bin mass-weighted VICE ages by [O/Fe]
    vice_grouped = group_by_bins(vice_stars, '[o/fe]', bins=ofe_bins)
    # weighted medians of VICE
    wm = lambda x: weighted_quantile(x, 'age', 'mstar', quantile=0.5)
    vice_medians = vice_grouped.apply(wm)
    # RMS of median difference
    notna = (pd.notna(apogee_medians) & pd.notna(vice_medians))
    median_diffs = vice_medians[notna] - apogee_medians[notna]
    d_rms = np.sqrt(np.average(median_diffs**2, weights=apogee_counts[notna]))
    
    return d_rms


def plot_vice_medians(ax, stars, ofe_lim=OFE_LIM, ofe_bin_width=OFE_BIN_WIDTH,
                      plot_low_mass_bins=True, low_mass_cutoff=0.01,
                      marker='s', small_marker='x', label='VICE', 
                      small_label=None, markersize=2):
    """
    Plot median stellar ages binned by [O/Fe] from VICE multizone data.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis on which to plot the medians.
    stars : pandas.DataFrame
        VICE multizone stars data, typically a subset of a galactic region.
    ofe_lim : tuple, optional
        Outermost bounds on [O/Fe]. The default is (-0.15, 0.55).
    ofe_bin_width : float, optional
        The [O/Fe] bin width in dex. The default is 0.05.
    plot_low_mass_bins : bool, optional
        Whether to include bins with little stellar mass. If True, these bins
        are shown with a different, smaller marker. The default is True.
    low_mass_cutoff : float, optional
        The cutoff as a fraction of the total mass for what constitutes a
        "low mass bin". The default is 0.01, meaning a bin which has less than
        1% the total mass of the simulation.
    marker : str, optional
        The marker style for the high-mass bins. The default is 's', a square.
    small_marker : str, optional
        The marker style for the low-mass bins. The default is 'x'.
    label : str, optional
        The main scatter plot / error bar label. The default is None.
    small_label : str, optional
        The small scatter plot / error bar label. The default is None.
    markersize : float, optional
        The marker size. The default is 2.
    """
    stars = stars.dropna(how='any')
    # Lambda functions for weighted quantiles
    wm = lambda x: weighted_quantile(x, 'age', 'mstar', quantile=0.5)
    wl = lambda x: weighted_quantile(x, 'age', 'mstar', quantile=0.16)
    wu = lambda x: weighted_quantile(x, 'age', 'mstar', quantile=0.84)
    # Define [O/Fe] bins
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    # Mass-weighted median and standard deviation of ages
    grouped = group_by_bins(stars, '[o/fe]', bins=ofe_bins)
    age_median = grouped.apply(wm)
    age_lower = grouped.apply(wl)
    age_upper = grouped.apply(wu)
    # Separate bins with very little mass and plot with different marker
    mtot = grouped['mstar'].sum()
    high_mass_bins = mtot[mtot >= low_mass_cutoff * mtot.sum()].index
    ax.errorbar(age_median[high_mass_bins], high_mass_bins, 
                xerr=(age_median[high_mass_bins] - age_lower[high_mass_bins], 
                      age_upper[high_mass_bins] - age_median[high_mass_bins]),
                yerr=ofe_bin_width/2,
                color='k', linestyle='none', capsize=1, elinewidth=0.5,
                capthick=0.5, marker=marker, markersize=markersize, label=label
    )
    # Plot bins with little (but non-zero) stellar mass with smaller markers
    if plot_low_mass_bins:
        low_mass_bins = mtot[(mtot > 0) & 
                             (mtot < low_mass_cutoff * mtot.sum())].index
        ax.errorbar(age_median[low_mass_bins], low_mass_bins, 
                    xerr=(age_median[low_mass_bins] - age_lower[low_mass_bins], 
                          age_upper[low_mass_bins] - age_median[low_mass_bins]),
                    yerr=ofe_bin_width/2,
                    color='k', linestyle='none', capsize=1, elinewidth=0.25,
                    capthick=0.25, marker=small_marker, markersize=markersize, 
                    markeredgewidth=0.5, label=small_label
        )
    return age_median, age_lower, age_upper, ofe_bins


def plot_astroNN_medians(ax, data, ofe_lim=OFE_LIM, ofe_bin_width=OFE_BIN_WIDTH,
                         plot_low_count_bins=True, low_count_cutoff=0.01,
                         marker='^', small_marker='2', label=None, 
                         small_label=None, markersize=2, age_col='ASTRONN_AGE'):
    """
    Plot median stellar ages binned by [O/Fe] from astroNN data.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis on which to plot the medians.
    data : pandas.DataFrame
        APOGEE data with astroNN ages, typically a subset of a galactic region.
    ofe_lim : tuple, optional
        Outermost bounds on [O/Fe]. The default is (-0.15, 0.55).
    ofe_bin_width : float, optional
        The [O/Fe] bin width in dex. The default is 0.05.
    plot_low_count_bins : bool, optional
        Whether to include bins with few stars. If True, these bins
        are shown with a different, smaller marker. The default is True.
    low_count_cutoff : float, optional
        The cutoff as a fraction of the total number of stars for what 
        constitutes a "low mass bin". The default is 0.01, meaning a bin which 
        has less than 1% the total number of stars.
    marker : str, optional
        The marker style for the high-count bins. The default is '^', an upwards
        pointing triangle.
    small_marker : str, optional
        The marker style for the low-count bins. The default is '2', an upwards
        pointing skinny triangle.
    label : str, optional
        The main scatter plot / error bar label. The default is None.
    small_label : str, optional
        The small scatter plot / error bar label. The default is None.
    markersize : float, optional
        The marker size. The default is 2.
    age_col : str, optional
        Name of column containing ages. The default is 'ASTRONN_AGE'
    """
    # Remove entries with no age estimate
    data.dropna(subset=age_col, inplace=True)
    # Define [O/Fe] bins
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    age_grouped = group_by_bins(data, 'O_FE', bins=ofe_bins)[age_col]
    age_median = age_grouped.median()
    age_lower = age_grouped.quantile(0.16)
    age_upper = age_grouped.quantile(0.84)
    # Separate bins with very little mass and plot with different markers
    counts = age_grouped.count()
    high_count_bins = counts[counts >= low_count_cutoff * counts.sum()].index
    ax.errorbar(age_median[high_count_bins], high_count_bins, 
                xerr=(age_median[high_count_bins] - age_lower[high_count_bins], 
                      age_upper[high_count_bins] - age_median[high_count_bins]),
                yerr=ofe_bin_width/2,
                color='r', linestyle='none', capsize=1, elinewidth=0.5,
                capthick=0.5, marker=marker, markersize=markersize, label=label,
    )
    # Plot bins with few (but non-zero) stars with a different, smaller marker
    if plot_low_count_bins:
        low_count_bins = counts[(counts > 0) & 
                                (counts < low_count_cutoff * counts.sum())].index
        ax.errorbar(age_median[low_count_bins], low_count_bins, 
                    xerr=(age_median[low_count_bins] - age_lower[low_count_bins], 
                          age_upper[low_count_bins] - age_median[low_count_bins]),
                    yerr=ofe_bin_width/2,
                    color='r', linestyle='none', capsize=1, elinewidth=0.25,
                    capthick=0.25, marker=small_marker, markersize=markersize, 
                    markeredgewidth=0.5, label=small_label,
        )
    return age_median, age_lower, age_upper, ofe_bins
        

def plot_feuillet2019(ax, galr_lim, absz_lim, **kwargs):
    r"""
    Plot the age-[X/Y] relation as reported by Feuillet et al. (2019) [1]_.

    Parameters
    ----------
    ax : ``axes``
        The matplotlib subplot to plot on.
    element_x : ``str`` [case-insensitive]
        The element X in age-[X/Y] relation.
    element_y : ``str`` [case-insensitive]
        The element Y in age-[X/Y] relation.
    min_rgal : ``float``
        Minimum galactocentric radius in kpc defining the region.
    max_rgal : ``float``
        Maximum galactocentric radius in kpc defining the region.
    min_absz : ``float``
        Minimum height above/below the disk midplane |z| in kpc defining the
        region.
    max_absz : ``float``
        Maximum height above/below the disk midplane |z| in kpc defining the
        region.
    label : ``bool`` [default : False]
        Whether or not to produce a legend handle for the plotted points with
        error bars.
    kwargs : varying types
        Additional keyword arguments to pass to ``pyplot.errorbar``.

    .. [1] Feuillet et al. (2019), MNRAS, 489, 1742
    """
    fname = 'ELEM_GAUSS_AGE_%02d_%02d_%02d_%02d_alpha.fits' % (
        galr_lim[0], galr_lim[1], 10 * absz_lim[0], 10 * absz_lim[1])
    datapath = paths.data / 'feuillet2019' / 'age_alpha' / fname
    age, ofe, age_disp, ofe_disp = feuillet2019_data(datapath)
    ax.errorbar(age, ofe, xerr=age_disp, yerr=ofe_disp,
                color='r', linestyle='none', capsize=1, elinewidth=0.5,
                capthick=0.5, marker='^', markersize=2, 
                label=AGE_LABELS['F19'], **kwargs)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='age_ofe.py',
        description='Generate a grid of [O/Fe] vs age plots comparing the' + \
            ' output of a VICE multizone run to APOGEE and astroNN data.'
    )
    parser.add_argument('evolution', metavar='EVOL',
                        help='Name of star formation history model')
    parser.add_argument('RIa', metavar='DTD',
                        help='Name of delay time distribution model')
    parser.add_argument('--migration', '-m', metavar='MIGR', 
                        choices=['diffusion', 'post-process', 'gaussian'], 
                        default='gaussian',
                        help='Name of migration prescription')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
                        default='winter',
                        help='Name of colormap for color-coding VICE ' + \
                             'output (default: winter)')
    parser.add_argument('-l', '--log', action='store_true',
                        help='Plot age on a log scale')
    parser.add_argument('-a', '--ages', choices=AGE_SOURCES, default='L23',
                        help='Source for age data (options: F19, M19, L23)')
    parser.add_argument('-u', '--uncertainties', action='store_true',
                        help='Model APOGEE uncertainties in VICE output')
    parser.add_argument('-s', '--score', action='store_true',
                        help='Numerically score the match between VICE ' + \
                             'output and data.')
    args = parser.parse_args()
    main(**vars(args))
