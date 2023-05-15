"""
This script plots a grid of panels showing [O/Fe] vs age in multiple galactic
regions. The plot includes a sample of VICE stellar populations, mass-weighted 
medians of the VICE output, and the median estimated ages from astroNN.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import vice
from scatter_plot_grid import setup_axes, setup_colorbar, plot_vice_sample
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    weighted_quantile, import_apogee, apogee_region, \
    group_by_bins, model_uncertainty, feuillet2019_data
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths

GALR_BINS = GALR_BINS[:-1]
AGE_LIM_LINEAR = (-1, 14)
AGE_LIM_LOG = (0.2, 20)
OFE_LIM = (-0.15, 0.55)
OFE_BIN_WIDTH = 0.05
AGE_SOURCES = ['F19', # Feuillet et al. 2018, solar neighborhood only
               'M19', # Mackereth et al. 2019, astroNN
               'L23'] # Leung et al. 2023, variational encoder-decoder
AGE_LABELS = {'F19': 'Feuillet et al. 2018',
              'M19': 'Mackereth et al. 2019',
              'L23': 'Leung et al. 2023'}

def main(evolution, RIa, migration='diffusion', verbose=False, cmap='winter',
         data_dir='../data/migration', log=False, ages='L23', 
         uncertainties=False):
    # Error handling
    if ages not in AGE_SOURCES:
        raise ValueError('Parameter "ages" must be in %s.' % AGE_SOURCES)
    # Import VICE multi-zone output data
    output_name = '/'.join(['diffusion', evolution, RIa])
    if verbose: 
        print('Importing VICE multizone data from %s/%s.vice' \
              % (data_dir, output_name))
    vice_stars = multioutput_to_pandas(output_name, data_dir)
    # Import APOGEE and astroNN data
    apogee_data = import_apogee(verbose=verbose)
    # Model uncertainties
    if uncertainties:
        # [O/Fe] uncertainty 
        ofe_err = apogee_data['O_FE_ERR'].median()
        vice_stars['[o/fe]'] = model_uncertainty(vice_stars['[o/fe]'], ofe_err)
        # age uncertainty
        if ages == 'F19':
            age_err = 0.15 # Feuillet et al. 2016 ApJ 817 40
            age_err_method = 'logarithmic'
        elif ages == 'M19':
            age_err = (apogee_data['ASTRONN_AGE_ERR'] / 
                       apogee_data['ASTRONN_AGE']).median()
            age_err_method = 'fractional'
        else: # L23
            age_err = apogee_data['LOG_LATENT_AGE_ERR'].median()
            age_err_method = 'logarithmic'
        vice_stars['age'] = model_uncertainty(vice_stars['age'], age_err,
                                              how=age_err_method)
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
    fig, axs = setup_axes(xlim=age_lim, ylim=OFE_LIM, 
                          xlabel='Age [Gyr]', ylabel='[O/Fe]',
                          galr_bins=local_galr_bins, absz_bins=local_absz_bins)
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]')
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
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
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, 
                                                   absz_lim, ZONE_WIDTH)
            plot_vice_sample(ax, vice_subset, 'age', '[o/fe]', 
                             cmap=cmap, norm=cbar.norm)
            # Plot Feuillet+ 2019 ages
            if ages == 'F19':
                plot_feuillet2019(ax, galr_lim, absz_lim)
            else:
                apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
                plot_astroNN_medians(ax, apogee_subset, 
                                     age_col=age_col, label=ages)
            plot_vice_medians(ax, vice_subset, label='VICE')
            # Add legend to top-right panel
            if i==0 and j==len(row)-1:
                ax.legend(loc='upper left', frameon=False, handletextpad=0.2)
                             
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
    fname = '%s_%s_%s.png' % (evolution, RIa, ages)
    fig.savefig(paths.figures / 'age_ofe' / fname, dpi=300)


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
    # Lambda functions for weighted quantiles
    wm = lambda x: weighted_quantile(x, 'age', 'mass', quantile=0.5)
    wl = lambda x: weighted_quantile(x, 'age', 'mass', quantile=0.16)
    wu = lambda x: weighted_quantile(x, 'age', 'mass', quantile=0.84)
    # Define [O/Fe] bins
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    # stars['odf_bin'] = pd.cut(stars['[o/fe]'], ofe_bins, 
    #                           labels=get_bin_centers(ofe_bins))
    # stars['mass_weighted_age'] = stars['age'] * stars['mass']
    # Mass-weighted median and standard deviation of ages
    # grouped = stars.groupby('odf_bin')
    grouped = group_by_bins(stars, '[o/fe]', bins=ofe_bins)
    age_median = grouped.apply(wm)
    age_lower = grouped.apply(wl)
    age_upper = grouped.apply(wu)
    # Separate bins with very little mass and plot with different marker
    mtot = grouped['mass'].sum()
    high_mass_bins = mtot[mtot >= low_mass_cutoff * mtot.sum()].index
    ax.errorbar(age_median[high_mass_bins], high_mass_bins, 
                xerr=(age_median[high_mass_bins] - age_lower[high_mass_bins], 
                      age_upper[high_mass_bins] - age_median[high_mass_bins]),
                yerr=ofe_bin_width/2,
                color='k', linestyle='none', capsize=1, elinewidth=0.5,
                capthick=0.5, marker=marker, markersize=markersize, label=label
    )
    # Plot bins with little stellar mass with smaller markers
    if plot_low_mass_bins:
        low_mass_bins = mtot[mtot < low_mass_cutoff * mtot.sum()].index
        ax.errorbar(age_median[low_mass_bins], low_mass_bins, 
                    xerr=(age_median[low_mass_bins] - age_lower[low_mass_bins], 
                          age_upper[low_mass_bins] - age_median[low_mass_bins]),
                    yerr=ofe_bin_width/2,
                    color='k', linestyle='none', capsize=1, elinewidth=0.25,
                    capthick=0.25, marker=small_marker, markersize=markersize, 
                    markeredgewidth=0.5, label=small_label
        )


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
    # Define [O/Fe] bins
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    # data['OFE_BIN'] = pd.cut(data['O_FE'], ofe_bins, 
    #                          labels=get_bin_centers(ofe_bins))
    # age_grouped = data.groupby('OFE_BIN')[age_col]
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
    # Plot bins with few stars with a different, smaller marker
    if plot_low_count_bins:
        low_count_bins = counts[counts < low_count_cutoff * counts.sum()].index
        ax.errorbar(age_median[low_count_bins], low_count_bins, 
                    xerr=(age_median[low_count_bins] - age_lower[low_count_bins], 
                          age_upper[low_count_bins] - age_median[low_count_bins]),
                    yerr=ofe_bin_width/2,
                    color='r', linestyle='none', capsize=1, elinewidth=0.25,
                    capthick=0.25, marker=small_marker, markersize=markersize, 
                    markeredgewidth=0.5, label=small_label,
        )
        

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
                label='F19', **kwargs)
    

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
                        choices=['diffusion', 'post-process'], 
                        default='diffusion',
                        help='Name of migration prescription')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
        default='winter',
        help='Name of colormap for color-coding VICE output (default: winter)')
    parser.add_argument('-l', '--log', action='store_true',
                        help='Plot age on a log scale')
    parser.add_argument('-a', '--ages', choices=AGE_SOURCES, default='L23',
                        help='Source for age data (options: F19, M19, L23)')
    parser.add_argument('-u', '--uncertainties', action='store_true',
                        help='Model APOGEE uncertainties in VICE output')
    args = parser.parse_args()
    main(**vars(args))
