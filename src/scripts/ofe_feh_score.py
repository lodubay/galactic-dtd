"""
This script scores how well a particular VICE model performes against APOGEE
data in the [O/Fe]-[Fe/H] plane. The score is based on the fraction of model
points which fall within the APOGEE 2-sigma contours. Galactic regions are
defined as they are in ofe_feh_vice.py and each region is given equal weight.
"""

import numpy as np
import matplotlib.pyplot as plt
from utils import import_allStar, multioutput_to_pandas, filter_multioutput_stars
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
from ofe_feh_vice import setup_axes, FEH_LIM, OFE_LIM
from ofe_feh_apogee import plot_contours, kde_path, read_kde
import paths

GALR_BINS = GALR_BINS[:-1]

def main(evolution, RIa, migration='diffusion', data_dir='../data/migration', 
         verbose=True, sigma=2):
    # Import VICE multi-zone output data
    output_name = '/'.join([migration, evolution, RIa])
    if verbose: 
        print('Importing VICE multizone data from %s/%s.vice' \
              % (data_dir, output_name))
    stars = multioutput_to_pandas(output_name, data_dir)
    
    # Test plot
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1)
    
    enclosed_fraction = []
    subset_mass = []
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            if verbose:
                print('RGal =', galr_lim, '|z| =', absz_lim)
            # Filter VICE stars by galactic region
            subset = filter_multioutput_stars(stars, galr_lim, absz_lim, 
                                              ZONE_WIDTH)
            # Import pre-generated 2D KDE of APOGEE data
            path = kde_path(galr_lim, absz_lim, savedir='../data/APOGEE/kde/ofe_feh/')
            xx, yy, logz = read_kde(path)
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
            vice_hist = np.histogram2d(subset['[fe/h]'], subset['[o/fe]'],
                                       bins=[xbins, ybins], density=True,
                                       weights=subset['mass'])[0]
            # Account for VICE stars which fall outside histogram range
            inrange = subset[(subset['[fe/h]'] >= xbins[0]) & 
                             (subset['[fe/h]'] < xbins[-1]) & 
                             (subset['[o/fe]'] >= ybins[0]) & 
                             (subset['[o/fe]'] < ybins[-1])]
            vice_hist *= inrange['mass'].sum() / subset['mass'].sum()
            axs[i,j].hist2d(xx.flatten(), yy.flatten(), bins=[xbins, ybins],
                           weights=vice_hist.flatten(), cmap='Greys')
            # fraction of VICE stars enclosed within APOGEE 2-sigma contour
            contour_frac = np.sum(vice_hist[kde_mask] * dx * dy)
            print(contour_frac)
            enclosed_fraction.append(contour_frac)
            subset_mass.append(np.sum(subset['mass']))
            # test plot
            vice_hist[~kde_mask] = 0.
            axs[i,j].hist2d(xx.flatten(), yy.flatten(), bins=[xbins, ybins],
                            weights=vice_hist.flatten(), cmap='Reds', cmin=1e-8)
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
    
    enclosed_fraction = np.array(enclosed_fraction)
    subset_mass = np.array(subset_mass)
    print('Mean unweighted enclosed fraction:\n', np.mean(enclosed_fraction))
    print('Mass-weighted mean enclosed fraction:\n', 
          np.average(enclosed_fraction, weights=subset_mass))
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    fig.savefig(paths.figures / 'ofe_feh_score.png', dpi=300)
    

if __name__ == '__main__':
    main('conroy22', 'plateau_width1000_slope11')
