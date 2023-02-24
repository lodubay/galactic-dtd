"""
This script scores how well different VICE models performe against APOGEE
data in the [O/Fe]-[Fe/H] plane. The score is based on the fraction of model
points which fall within the APOGEE 2-sigma contours. Galactic regions are
defined as they are in ofe_feh_vice.py and each region is given equal weight.
An alternative score weights each region by the APOGEE sample size.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import import_apogee, multioutput_to_pandas, \
    filter_multioutput_stars, apogee_region
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
from ofe_feh_vice import setup_axes, FEH_LIM, OFE_LIM
from ofe_feh_apogee import gen_kde
import paths

GALR_BINS = GALR_BINS[:-1]
SIGMA = 2

SFHs = ['insideout', 'lateburst', 'conroy22', 'twoinfall']
DTDs = ['powerlaw_slope11', 
        'powerlaw_slope14', 
        'exponential_timescale15',
        'exponential_timescale30',
        'plateau_width300_slope11',
        'plateau_width1000_slope11',
        'prompt_peak050_stdev015_timescale30',
        'triple_delay040']

def main(verbose=True, sigma=SIGMA):
    if verbose:
        print('Importing APOGEE allStar data...')
    apogee_data = import_apogee()
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['Unweighted', 'Weighted'],
    )
    
    for evolution in SFHs:
        for RIa in DTDs:
            output_name = '/'.join(('diffusion', evolution, RIa))
            scores = score_output(output_name, apogee_data, sigma=sigma, verbose=verbose)
            summary_table.loc[RIa, evolution] = np.array(scores).round(3)
    
    if verbose: print(summary_table)
    summary_table.to_csv(paths.figures / 'ofe_feh/scores/scores.csv')


def score_output(output_name, apogee_data, sigma=SIGMA, verbose=True):
    """
    Score the given VICE multizone output in the [O/Fe]-[Fe/H] plane by the 
    fraction of stellar populations enclosed within the APOGEE 2-sigma contour.
    
    Parameters
    ----------
    output_name : str
        VICE multizone output name, typically '<migration>/<evolution>/<RIa>'
    apogee_data : pandas.DataFrame
        APOGEE data
    sigma : float, optional
        Standard deviation APOGEE contour to test against VICE. The default
        is 2.
    verbose : bool, optional
        If True, print diagnostic info to terminal.
    
    Returns
    -------
    unweighted : float
        Mean unweighted enclosed fraction of stars over all Galactic regions
    weighted : float
        Mean enclosed fraction, weighted by total mass of VICE stars in
        each region
    """
    
    # Import VICE multi-zone output data
    if verbose: 
        print('Importing VICE multizone data from %s/%s.vice' \
              % ((str(paths.migration)), output_name))
    stars = multioutput_to_pandas(output_name)
    
    # Set up diagnostic plot
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1)
    
    enclosed_fraction = []
    weights = []
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            # Filter VICE stars by galactic region
            vice_subset = filter_multioutput_stars(stars, galr_lim, absz_lim, 
                                              ZONE_WIDTH)
            # Filter APOGEE data by galactic region
            apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
            # Generate or import 2D KDE of APOGEE data
            xx, yy, logz = gen_kde(apogee_data, absz_lim=absz_lim, 
                                   galr_lim=galr_lim)
            scaled_density = np.exp(logz) / np.max(np.exp(logz))
            # density threshold for N-sigma contour
            density_threshold = np.exp(-0.5 * sigma**2)
            kde_mask = (scaled_density > density_threshold)
            # Convert 2D KDE coordinates to bins in x and y
            dx = xx[1,0] - xx[0,0]
            xbins = np.arange(xx[0,0]-dx/2, xx[-1,0]+dx, dx)
            dy = yy[0,1] - yy[0,0]
            ybins = np.arange(yy[0,0]-dy/2, yy[0,-1]+dy, dy)
            # Bin VICE stars according to APOGEE 2D KDE histogram
            vice_hist = np.histogram2d(vice_subset['[fe/h]'], 
                                       vice_subset['[o/fe]'],
                                       bins=[xbins, ybins], density=True,
                                       weights=vice_subset['mass'])[0]
            # Account for VICE stars which fall outside histogram range
            inrange = vice_subset[(vice_subset['[fe/h]'] >= xbins[0]) & 
                                  (vice_subset['[fe/h]'] < xbins[-1]) & 
                                  (vice_subset['[o/fe]'] >= ybins[0]) & 
                                  (vice_subset['[o/fe]'] < ybins[-1])]
            vice_hist *= inrange['mass'].sum() / vice_subset['mass'].sum()
            vice_hist_max = vice_hist.max()
            axs[i,j].hist2d(xx.flatten(), yy.flatten(), bins=[xbins, ybins],
                            weights=vice_hist.flatten(), cmap='Greys')
            # fraction of VICE stars enclosed within APOGEE 2-sigma contour
            contour_frac = np.sum(vice_hist[kde_mask] * dx * dy)
            enclosed_fraction.append(contour_frac)
            # Weight by number of APOGEE stars in region
            weights.append(apogee_subset.shape[0])
            # test plot
            vice_hist[~kde_mask] = 0.
            axs[i,j].hist2d(xx.flatten(), yy.flatten(), bins=[xbins, ybins],
                            weights=vice_hist.flatten(), cmap='Reds', 
                            cmin=1e-8, vmax=vice_hist_max)
            axs[i,j].contour(xx, yy, scaled_density, [density_threshold],
                             linewidths=0.5)
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.1, 0.1, r'$%s\leq |z| < %s$ kpc' % absz_lim,
                        transform=ax.transAxes, size=8)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
            # Label enclosed fraction
            ax.text(0.9, 0.9, r'$f_{\rm{enc}}=%s$' % round(contour_frac, 2),
                    transform=ax.transAxes, size=8, ha='right', va='top')
    
    unweighted = np.mean(enclosed_fraction)
    if verbose:
        print('Mean unweighted enclosed fraction:\n', unweighted)
    weighted = np.average(enclosed_fraction, weights=weights)
    if verbose:
        print('Mass-weighted mean enclosed fraction:\n', weighted)
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    figname = '_'.join(output_name.split('/')[1:])
    fig.savefig(paths.figures / ('ofe_feh/scores/%s.png' % figname), dpi=300)
    plt.close()
    
    del stars
    return unweighted, weighted
    

if __name__ == '__main__':
    main()
