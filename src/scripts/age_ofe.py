"""
Plot a grid of panels showing [O/Fe] vs age in various galactic regions.
Includes a sample of VICE stellar populations, mass-weighted medians of the
VICE output, and the median observed values from APOGEE and astroNN.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice
from scatter_plot_grid import setup_axes, setup_colorbar
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    sample_dataframe, get_bin_centers, weighted_quantile, import_astroNN, \
    apogee_region
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths


GALR_BINS = GALR_BINS[:-1]
AGE_LIM = (0.2, 20)
OFE_LIM = (-0.15, 0.55)


def main(evolution, RIa, migration='diffusion', verbose=False, cmap='winter',
         data_dir='../data/migration'):
    # Import VICE multi-zone output data
    vice_stars = multioutput_to_pandas('/'.join(['diffusion', evolution, RIa]),
                                       data_dir)
    
    # Set up figure
    fig, axs = setup_axes(xlim=AGE_LIM, ylim=OFE_LIM, 
                          xlabel='Age [Gyr]', ylabel='[O/Fe]')
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]')
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    # Plot sampled points and medians
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, 
                                                   absz_lim, ZONE_WIDTH)
            plot_vice_sample(ax, vice_subset, cmap=cmap, norm=cbar.norm)
    
    # Output figure
    fname = '%s_%s.png' % (evolution, RIa)
    fig.savefig(paths.figures / 'age_ofe' / fname, dpi=300)


def plot_vice_sample(ax, stars, cmap=None, norm=None, sampled=True, 
                     nsamples=10000, markersize=0.1):
    """
    Generate a scatter plot of [O/Fe] vs age from VICE multizone stars data.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes on which to draw the scatter plot.
    stars : pandas.DataFrame
        VICE multizone stars data, usually a subset defined by a particular
        Galactic region.
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap to color-code the scattered points. The default is None.
    norm : matplotlib.colors.Normalize, optional
        Normalization of the color-bar mapping. The default is None.
    sampled : bool, optional
        If True, sample the VICE stars (weighted by mass) instead of plotting 
        them all. The default is True.
    nsamples : int, optional
        Number of stellar populations to sample from the VICE output. The
        default is 10000.
    markersize : float, optional
        Size of scatter plot markers. The default is 0.1.
    """
    if sampled:
        # weight random sample by particle mass
        sample_weights = stars['mass'] / stars['mass'].sum()
        sample = sample_dataframe(stars, nsamples, weights=sample_weights)
    else:
        sample = stars.copy()
    # Scatter plot of stellar particles
    ax.scatter(sample['age'], sample['[o/fe]'], s=markersize,
               c=sample['zone_origin'] * ZONE_WIDTH, cmap=cmap,
               norm=norm, rasterized=True, edgecolor='none')


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
    args = parser.parse_args()
    main(**vars(args))
