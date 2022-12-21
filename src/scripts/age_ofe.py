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
    get_bin_centers, weighted_quantile, import_astroNN, apogee_region
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths

GALR_BINS = GALR_BINS[:-1]
AGE_LIM_LINEAR = (-1, 14)
AGE_LIM_LOG = (0.2, 20)
OFE_LIM = (-0.15, 0.55)
OFE_BIN_WIDTH = 0.05

def main(evolution, RIa, migration='diffusion', verbose=False, cmap='winter',
         data_dir='../data/migration', log=False):
    # Import VICE multi-zone output data
    output_name = '/'.join(['diffusion', evolution, RIa])
    if verbose: 
        print('Importing VICE multizone data from %s/%s.vice' \
              % (data_dir, output_name))
    vice_stars = multioutput_to_pandas(output_name, data_dir)
    # Import APOGEE and astroNN data
    astroNN_data = import_astroNN(verbose=verbose)
    
    # Set x-axis limits
    if log:
        age_lim = AGE_LIM_LOG
    else:
        age_lim = AGE_LIM_LINEAR
    
    # Set up figure
    fig, axs = setup_axes(xlim=age_lim, ylim=OFE_LIM, 
                          xlabel='Age [Gyr]', ylabel='[O/Fe]')
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]')
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    # Plot sampled points and medians
    if verbose: 
        print('Plotting [O/Fe] vs age in galactic regions...')
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            if verbose:
                print('\tRGal=%s kpc, |z|=%s kpc' \
                      % (str(galr_lim), str(absz_lim)))
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, 
                                                   absz_lim, ZONE_WIDTH)
            plot_vice_sample(ax, vice_subset, 'age', '[o/fe]', 
                             cmap=cmap, norm=cbar.norm)
            plot_vice_medians(ax, vice_subset.copy())
            astroNN_subset = apogee_region(astroNN_data, galr_lim, absz_lim)
            plot_astroNN_medians(ax, astroNN_subset.copy())
                             
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
    fname = '%s_%s.png' % (evolution, RIa)
    fig.savefig(paths.figures / 'age_ofe' / fname, dpi=300)


def plot_vice_medians(ax, stars, ofe_lim=OFE_LIM, ofe_bin_width=OFE_BIN_WIDTH,
                      plot_low_mass_bins=True, low_mass_cutoff=0.01,
                      marker='s', small_marker='x', label=None, 
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
    stars['odf_bin'] = pd.cut(stars['[o/fe]'], ofe_bins, 
                              labels=get_bin_centers(ofe_bins))
    stars['mass_weighted_age'] = stars['age'] * stars['mass']
    # Mass-weighted median and standard deviation of ages
    grouped = stars.groupby('odf_bin')
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
                         small_label=None, markersize=2):
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
    """
    # Define [O/Fe] bins
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    data['OFE_BIN'] = pd.cut(data['O_FE'], ofe_bins, 
                             labels=get_bin_centers(ofe_bins))
    age_grouped = data.groupby('OFE_BIN')['ASTRONN_AGE']
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
    args = parser.parse_args()
    main(**vars(args))
