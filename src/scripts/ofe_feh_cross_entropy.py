"""
Calculate the cross entropy in the [O/Fe]-[Fe/H] plane between VICE stellar
populations from multizone runs and APOGEE.
"""

import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import entropy
from ofe_feh_apogee import gen_kde as apogee_kde
from utils import cross_entropy, kde2D, filter_multioutput_stars, \
    apogee_region, multioutput_to_pandas, import_apogee, sample_dataframe
from ofe_feh_vice import setup_axes, FEH_LIM, OFE_LIM
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths

SFHs = ['insideout', 'lateburst', 'conroy22', 'twoinfall']
DTDs = ['powerlaw_slope11', 
        'powerlaw_slope14', 
        'exponential_timescale15',
        'exponential_timescale30',
        'plateau_width300_slope11',
        'plateau_width1000_slope11',
        'prompt_peak050_stdev015_timescale30',
        'triple_delay040']

def main(verbose=True):
    output_name = 'diffusion/insideout/powerlaw_slope11'
    apogee_data = import_apogee(verbose=verbose)
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['Unweighted', 'Weighted'],
    )
    
    for evolution in SFHs:
        for RIa in DTDs:
            output_name = '/'.join(('diffusion', evolution, RIa))
            scores = ce_regions(output_name, apogee_data, verbose=verbose)
            summary_table.loc[RIa, evolution] = np.array(scores).round(2)
    
    if verbose: print(summary_table)
    summary_table.to_csv(paths.figures / 'ofe_feh/cross_entropy/scores.csv')


def ce_regions(output_name, apogee_data, galr_bins=GALR_BINS[:-1],
               absz_bins=ABSZ_BINS, verbose=False, bandwidth=0.02):
    """
    Calculate the cross entropy for each galactic region.
    
    Parameters
    ----------
    output_name : str
        VICE multizone output name, typically '<migration>/<evolution>/<RIa>'
    apogee_data : pandas.DataFrame
        APOGEE data containing relevant columns.
    verbose : bool, optional
        If True, print diagnostic info to terminal.
    bandwidth : float, optional
        Bandwidth of the KDE. The default is 0.02.

    Returns
    -------
    cedf : pandas.DataFrame
        DataFrame containing cross entropy in each region. Columns are labeled
        by a tuple of Rgal bounds, and rows by a tuple of |z| bounds.
    unweighted : float
        Mean unweighted enclosed fraction of stars over all Galactic regions
    weighted : float
        Mean enclosed fraction, weighted by total number of APOGEE stars
        in each region
    """
    vice_stars = multioutput_to_pandas(output_name, verbose=verbose)
    # Initialize figure
    fig, axs = setup_axes(len(absz_bins)-1, len(galr_bins)-1)
    
    # initialize dataframe
    cols = [(galr_bins[i], galr_bins[i+1]) for i in range(len(galr_bins)-1)]
    rows = [(absz_bins[i], absz_bins[i+1]) for i in range(len(absz_bins)-1)]
    rows.reverse()
    cedf = pd.DataFrame([], index=rows, columns=cols)
    
    # Iterate through galactic regions
    weights = []
    ce_list = []
    for i, absz_lim in enumerate(rows):
        for j, galr_lim in enumerate(cols):
            if verbose:
                print('|z| =', absz_lim, ', Rgal =', galr_lim)
            # Filter APOGEE data by galactic region
            apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
            # Weight by number of APOGEE stars in region
            weights.append(apogee_subset.shape[0])
            # Import or generate 2D KDE of APOGEE data
            xx, yy, apogee_logz = apogee_kde(apogee_subset, bandwidth=bandwidth,
                                             absz_lim=absz_lim, galr_lim=galr_lim)
            # Filter VICE stars by galactic region
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, 
                                                   absz_lim, ZONE_WIDTH)
            # KDE of VICE subset
            xx, yy, vice_logz = kde2D(vice_subset['[fe/h]'], vice_subset['[o/fe]'],
                                      bandwidth, xbins=xx, ybins=yy)
            ce = cross_entropy(np.exp(apogee_logz), np.exp(vice_logz))
            cedf.iloc[i, j] = ce
            ce_list.append(ce)
            if verbose:
                print('\tCE =', ce)
            
            # Plot
            ax = axs[i,j]
            # APOGEE contour levels
            levels = -0.5 * np.array([2, 1])**2
            linestyles = ['--', '-']
            ax.contour(xx, yy, apogee_logz - np.max(apogee_logz), levels, 
                       colors='r', linestyles=linestyles, linewidths=0.5)
            # VICE contour levels
            ax.contour(xx, yy, vice_logz - np.max(vice_logz), levels, 
                       colors='k', linestyles=linestyles, linewidths=0.5)
            # VICE scatterplot
            # weight random sample by particle mass
            sample_weights = vice_subset['mass'] / vice_subset['mass'].sum()
            sample = sample_dataframe(vice_subset, 10000, weights=sample_weights)
            # Scatter plot of random sample of stellar particles
            ax.scatter(sample['[fe/h]'], sample['[o/fe]'], s=0.1,
                       rasterized=True, edgecolor='none')
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.1, r'$%s\leq |z| < %s$ kpc' % absz_lim,
                        transform=ax.transAxes, size=8)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
            # Label cross entropy
            ax.text(0.9, 0.9, r'CE=%.02f' % ce,
                    transform=ax.transAxes, size=8, ha='right', va='top')
    
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    figname = '_'.join(output_name.split('/')[1:])
    fig.savefig(paths.figures / ('ofe_feh/cross_entropy/%s.png' % figname), 
                dpi=300)
    plt.close()
    
    unweighted = np.mean(ce_list)
    if verbose:
        print('Mean unweighted cross-entropy:\n', unweighted)
    weighted = np.average(ce_list, weights=weights)
    if verbose:
        print('Mass-weighted mean cross-entropy:\n', weighted)
    
    # return cedf
    return unweighted, weighted

if __name__ == '__main__':
    main()
