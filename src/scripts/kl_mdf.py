"""
This script calculates the KL divergence between the APOGEE MDF and a VICE
multizone model output.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from feh_distribution import apogee_mdf, vice_mdf
from utils import kl_divergence, multioutput_to_pandas, axes_grid, get_bin_centers
from apogee_tools import import_apogee, apogee_region
from _globals import GALR_BINS, ABSZ_BINS
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

FEH_LIM = (-2., 1.)
BIN_WIDTH = 0.01 # width of MDF bins
SMOOTHING = 0.2 # width of boxcar smoothing

def main(verbose=True):
    apogee_data = import_apogee(verbose=verbose)
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['Weighted KL Divergence'],
    )
    
    for evolution in SFHs:
        for RIa in DTDs:
            output_name = '/'.join(('diffusion', evolution, RIa))
            kld_mean = mean_kld(output_name, apogee_data, xlim=FEH_LIM, 
                                bin_width=BIN_WIDTH, smooth_width=SMOOTHING, 
                                verbose=verbose, diagnostic=True)
            summary_table.loc[RIa, evolution] = np.array(kld_mean).round(3)
    
    if verbose: print(summary_table)
    summary_table.to_csv(paths.figures / 'mdf_feh/kl_divergence/scores.csv')


def mean_kld(output_name, apogee_data, weighted=True, verbose=False,
             diagnostic=False, **kwargs):
    """
    Calculate the mean KL divergence in all regions.
    
    Parameters
    ----------
    output_name : str
        VICE multizone output name, typically '<migration>/<evolution>/<RIa>'
    apogee_data : pandas.DataFrame
        APOGEE data containing relevant columns.
    weighted : bool, optional
        If True, weight the mean by the number of APOGEE stars in each region.
        The default is True.
    verbose : bool, optional
        If True, print verbose output. The default is False.
    diagnostic : bool, optional
        If True, create a diagnostic plot which is saved in 
        figures/mdf_feh/kld/. The default is False.
        
    **kwargs passed to vice_mdf() and apogee_mdf()
    xlim : tuple, optional
        Minimum and maximum abundance values. The default is (-1.1, 0.6).
    bin_width : float, optional
        Width of histogram bins in x-axis units. The default is 0.01.
    smooth_width : float, optional
        Width of boxcar smoothing in x-axis units. The default is 0.2.
    
    Returns
    -------
    kl_mean : float
        Mean (weighted or unweighted) KL divergence.
    """
    vice_stars = multioutput_to_pandas(output_name, verbose=verbose)
    
    # Initialize figure
    if diagnostic:
        fig, axs = axes_grid(len(ABSZ_BINS)-1, len(GALR_BINS)-1, xlim=FEH_LIM, ylim=(0, 3))
    
    kld_list = []
    weights = []
    for i in range(len(ABSZ_BINS) - 1):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j in range(len(GALR_BINS) - 1):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            if verbose:
                print('|z| =', absz_lim, ', Rgal =', galr_lim)
            
            vice_dist, bin_edges = vice_mdf(vice_stars,'[fe/h]', galr_lim, 
                                            absz_lim, **kwargs)
            apogee_dist, bin_edges = apogee_mdf(apogee_data, 'FE_H', galr_lim, 
                                                absz_lim, **kwargs)
            dx = bin_edges[1] - bin_edges[0]
            kld = kl_divergence(apogee_dist, vice_dist, dx)
            kld_list.append(kld)
            if verbose:
                print('\tKL divergence =', kld)
            
            if weighted:
                apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
                apogee_subset.dropna(axis=0, subset=['FE_H'], inplace=True)
                weights.append(apogee_subset.shape[0])
            else:
                weights.append(1.)
                
            # Plot distributions
            if diagnostic:
                ax = axs[i,j]
                xarr = get_bin_centers(bin_edges)
                ax.plot(xarr, vice_dist, 'k-', label='VICE', linewidth=1)
                ax.plot(xarr, apogee_dist, 'r--', label='APOGEE', linewidth=1)
                # Label axes
                if i == len(ABSZ_BINS)-2:
                    ax.set_xlabel('[Fe/H]')
                if j == 0:
                    ax.set_ylabel('p([Fe/H])')
                    ax.text(0.1, 0.5, r'$%s\leq |z| < %s$ kpc' % absz_lim,
                            transform=ax.transAxes, size=8, va='center', ha='left', rotation=90)
                if i == 0:
                    ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
                # Label cross entropy
                ax.text(0.92, 0.92, r'KL=%.03f' % kld,
                        transform=ax.transAxes, size=8, ha='right', va='top')
                
    # Finish plot
    if diagnostic:
        axs[0,-1].legend(frameon=False, loc='upper left', handlelength=1.3)
        figname = '_'.join(output_name.split('/')[1:])
        plt.savefig(paths.figures / ('mdf_feh/kl_divergence/%s.png' % figname),
                    dpi=300)
        plt.close()
                
    kld_mean = np.average(kld_list, weights=weights)
    if verbose:
        print('Mean weighted KL divergence:\n', kld_mean, '\n')
        
    return kld_mean


if __name__ == '__main__':
    main()
