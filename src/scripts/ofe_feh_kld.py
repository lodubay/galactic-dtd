"""
Calculate the 2D Kullback-Leibler (KL) divergence in the [O/Fe]-[Fe/H] plane 
between VICE stellar populations from multizone runs and APOGEE.
"""

import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
from utils import kl_div_2D, filter_multioutput_stars, \
    apogee_region, multioutput_to_pandas, import_apogee, sample_dataframe, \
    model_uncertainties
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
    apogee_data = import_apogee(verbose=verbose)
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['Unweighted', 'Weighted'],
    )
    
    for evolution in SFHs:
        for RIa in DTDs:
            output_name = '/'.join(('diffusion', evolution, RIa))
            scores = kld_regions(output_name, apogee_data, verbose=verbose)
            summary_table.loc[RIa, evolution] = np.array(scores).round(2)
    
    if verbose: print(summary_table)
    summary_table.to_csv(paths.figures / 'ofe_feh/kld/scores.csv')


def kld_regions(output_name, apogee_data, galr_bins=GALR_BINS[:-1],
               absz_bins=ABSZ_BINS, verbose=False):
    """
    Calculate the KL divergence for each galactic region.
    
    Parameters
    ----------
    output_name : str
        VICE multizone output name, typically '<migration>/<evolution>/<RIa>'
    apogee_data : pandas.DataFrame
        APOGEE data containing relevant columns.
    verbose : bool, optional
        If True, print diagnostic info to terminal.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing KL divergence in each region. Columns are labeled
        by a tuple of Rgal bounds, and rows by a tuple of |z| bounds.
    unweighted : float
        Mean unweighted KL divergence of stars over all Galactic regions
    weighted : float
        Mean KL divergence, weighted by total number of APOGEE stars
        in each region
    """
    vice_stars = multioutput_to_pandas(output_name, verbose=verbose)
    vice_stars = model_uncertainties(vice_stars.copy())
    # Initialize figure
    fig, axs = setup_axes(len(absz_bins)-1, len(galr_bins)-1)
    
    # initialize dataframe
    cols = [(galr_bins[i], galr_bins[i+1]) for i in range(len(galr_bins)-1)]
    rows = [(absz_bins[i], absz_bins[i+1]) for i in range(len(absz_bins)-1)]
    rows.reverse()
    # cedf = pd.DataFrame([], index=rows, columns=cols)
    
    # Iterate through galactic regions
    weights = []
    kld_list = []
    for i, absz_lim in enumerate(rows):
        for j, galr_lim in enumerate(cols):
            if verbose:
                print('|z| =', absz_lim, ', Rgal =', galr_lim)
            # Filter APOGEE data by galactic region
            apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
            apogee_subset.dropna(subset=['FE_H', 'O_FE'], inplace=True)
            # Weight by number of APOGEE stars in region
            weights.append(apogee_subset.shape[0])
            # Filter VICE stars by galactic region
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, 
                                                   absz_lim, ZONE_WIDTH)
            # 2D KL divergence between APOGEE (true) and VICE (approximate)
            kld = kl_div_2D(apogee_subset[['FE_H', 'O_FE']], 
                            vice_subset[['[fe/h]', '[o/fe]']])
            # cedf.iloc[i, j] = ce
            kld_list.append(kld)
            if verbose:
                print('\tKL divergence =', kld)
            
            # Plot
            ax = axs[i,j]
            # random sample of APOGEE stars
            apogee_sample = sample_dataframe(apogee_subset, 2000)
            ax.scatter(apogee_sample['FE_H'], apogee_sample['O_FE'], 
                       s=0.1, c='r', rasterized=True, edgecolor='none',
                       label='APOGEE')
            # VICE scatterplot
            # weight random sample by particle mass
            sample_weights = vice_subset['mass'] / vice_subset['mass'].sum()
            vice_sample = sample_dataframe(vice_subset, 10000, 
                                           weights=sample_weights)
            # Scatter plot of random sample of stellar particles
            ax.scatter(vice_sample['[fe/h]'], vice_sample['[o/fe]'], 
                       s=0.1, c='k', rasterized=True, edgecolor='none',
                       label='VICE')
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.1, r'$%s\leq |z| < %s$ kpc' % absz_lim,
                        transform=ax.transAxes, size=8)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
            # Label KL divergence
            ax.text(0.9, 0.9, r'KLD=%.02f' % kld,
                    transform=ax.transAxes, size=8, ha='right', va='top')
            # Add legend to top-right panel
            if i==0 and j==len(cols)-1:
                ax.legend(loc='lower left', frameon=False,
                          borderpad=0., handletextpad=0.2, markerscale=5.)
    
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    figname = '_'.join(output_name.split('/')[1:])
    fig.savefig(paths.figures / ('ofe_feh/kld/%s.png' % figname), 
                dpi=300)
    plt.close()
    
    unweighted = np.mean(kld_list)
    if verbose:
        print('Mean unweighted 2D KL divergence:\n', unweighted)
    weighted = np.average(kld_list, weights=weights)
    if verbose:
        print('Mass-weighted mean 2D KL divergence:\n', weighted)
    
    # return cedf
    return unweighted, weighted

if __name__ == '__main__':
    main()
