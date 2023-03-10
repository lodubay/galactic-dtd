"""
Test the cross-entropy between two samples drawn from the same 2D distribution.
"""

import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import entropy
from utils import cross_entropy, kde2D, \
    apogee_region, import_apogee, sample_dataframe
from ofe_feh_vice import setup_axes, FEH_LIM, OFE_LIM
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths

def main(verbose=True):
    # Import APOGEE data
    data = import_apogee(verbose=verbose)
    
    # Initialize figure
    galr_bins = GALR_BINS[:-1]
    absz_bins = ABSZ_BINS
    fig, axs = setup_axes(len(absz_bins)-1, len(galr_bins)-1)
    
    # initialize cross-entropy dataframe
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
            subset = apogee_region(data, galr_lim, absz_lim)
            subset = subset.copy().dropna(axis=0, subset=['FE_H', 'O_FE'])
            # Weight by number of APOGEE stars in region
            weights.append(subset.shape[0])
            # Split into 2 random subsets
            sample1 = sample_dataframe(subset, int(subset.shape[0]/2), reset=False)
            sample2 = subset.drop(sample1.index, axis=0).reset_index(drop=True)
            sample1.reset_index(inplace=True, drop=True)
            # Generate KDEs
            bandwidth = 0.02
            xx, yy, logz1 = kde2D(sample1['FE_H'], sample1['O_FE'], bandwidth)
            xx, yy, logz2 = kde2D(sample2['FE_H'], sample2['O_FE'], bandwidth,
                                  xbins=xx, ybins=yy)
            dx = xx[1,0] - xx[0,0]
            dy = yy[0,1] - yy[0,0]
            
            ce = cross_entropy(np.exp(logz1), np.exp(logz2))
            cedf.iloc[i, j] = ce
            ce_list.append(ce)
            if verbose:
                print('\tCE =', ce)
            
            # Plot
            ax = axs[i,j]
            # APOGEE contour levels
            levels = -0.5 * np.array([2, 1])**2
            linestyles = ['--', '-']
            ax.contour(xx, yy, logz1 - np.max(logz1), levels,
                       colors='r', linestyles=linestyles, linewidths=0.5)
            # VICE contour levels
            ax.contour(xx, yy, logz2 - np.max(logz2), levels,
                       colors='k', linestyles=linestyles, linewidths=0.5)
            
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
    
    fig.savefig(paths.figures / 'cross_entropy_test.png', dpi=300)
    plt.close()
    
    unweighted = np.mean(ce_list)
    if verbose:
        print('Mean unweighted cross-entropy:\n', unweighted)
    weighted = np.average(ce_list, weights=weights)
    if verbose:
        print('Mass-weighted mean cross-entropy:\n', weighted)

if __name__ == '__main__':
    main()
