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
from ofe_feh_apogee import plot_contours, kde_path, read_kde

def main(evolution, RIa, migration='diffusion', data_dir='../data/migration', 
         verbose=True, sigma=2):
    # Import VICE multi-zone output data
    output_name = '/'.join([migration, evolution, RIa])
    if verbose: 
        print('Importing VICE multizone data from %s/%s.vice' \
              % (data_dir, output_name))
    stars = multioutput_to_pandas(output_name, data_dir)
    
    total_enclosed_fraction = 0.
    for i in range(len(ABSZ_BINS)-1):
        absz_lim = (ABSZ_BINS[i], ABSZ_BINS[i+1])
        for j in range(len(GALR_BINS)-2):
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
            vice_hist = np.histogram2d(stars['[fe/h]'], stars['[o/fe]'],
                                       bins=[xbins, ybins], density=True)[0]
            # fraction of VICE stars enclosed within APOGEE 2-sigma contour
            contour_frac = np.sum(vice_hist[kde_mask] * dx * dy)
            print(contour_frac)
            total_enclosed_fraction += contour_frac
    
    nregions = (len(ABSZ_BINS)-1) * (len(GALR_BINS)-1)
    mean_enclosed_fraction = total_enclosed_fraction / nregions
    print(mean_enclosed_fraction)
    

if __name__ == '__main__':
    main('insideout', 'powerlaw_slope14')
