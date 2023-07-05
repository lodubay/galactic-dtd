"""
This script plots the distributions of migration distance as a function of
stellar population age and formation radius.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from multizone_stars import MultizoneStars
from utils import box_smooth, get_color_list
import paths
from _globals import TWO_COLUMN_WIDTH

MIGRATION_SCHEMES = ['diffusion', 'gaussian']
ROW_LABELS = ['h277 analog', 'Gaussian sampling']
RFORM_BINS = [3, 5, 7, 9, 11, 13]
AGE_BINS = [0, 2, 4, 6, 8, 10, 12]

def main():
    rfinal_bins = np.arange(0, 20.01, 0.01)
    
    # pick discrete colors
    cmap = plt.get_cmap('jet')
    colors = get_color_list(cmap, AGE_BINS)
    
    fig, axs = plt.subplots(len(MIGRATION_SCHEMES), len(RFORM_BINS), 
                            figsize=(TWO_COLUMN_WIDTH, 3), 
                            sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0., hspace=0.)
    
    for m, row in enumerate(axs):
        mig = MIGRATION_SCHEMES[m]
        output_name = '/'.join((mig, 'insideout/powerlaw_slope11/diskmodel'))
        mzs = MultizoneStars.from_output(output_name)
        row[2].text(0.9, 0.9, ROW_LABELS[m], 
                    ha='right', transform=row[2].transAxes)
        for i, ax in enumerate(row):
            rform_lim = tuple(RFORM_BINS[i:i+2])
            ax.vlines(rform_lim, ymin=0, ymax=0.6, linestyle=':', color='k')
            if m == 0:
                ax.set_title(r'$%d\leq R_{\rm{form}} < %d$ kpc' % rform_lim)
            for j, color in enumerate(colors):
                age_lim = tuple(AGE_BINS[j:j+2])
                # bin data by formation radius and age
                subset = mzs.filter({'age': age_lim, 'galr_origin': rform_lim})
                # limit to bins with a meaningful number of stars
                if subset.stars.shape[0] > 100:
                    hist, _ = np.histogram(subset('galr_final'), rfinal_bins, 
                                           density=True, weights=subset('mstar'))
                    # apply boxcar smoothing with a width of 0.5 kpc
                    hist_smooth = box_smooth(hist, rfinal_bins, 0.5)
                    bin_centers = (rfinal_bins[:-1] + rfinal_bins[1:]) / 2
                    ax.plot(bin_centers, hist_smooth, c=color, ls='-', 
                            label=r'$%d - %d$ Gyr' % age_lim)
    
    axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    axs[0,0].set_xlim((-2, 22))
    axs[0,0].set_ylim((0, 0.45))
    for ax in axs[-1]:
        ax.set_xlabel(r'$R_{\rm{final}}$ [kpc]')
    for ax in axs[:,0]:
        ax.set_ylabel('PDF')
    axs[0,0].legend(title='Age bins', loc='upper right', fontsize=8, frameon=False)
    plt.savefig(paths.figures / 'radial_migration.png', dpi=300)

if __name__ == '__main__':
    main()
