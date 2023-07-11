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
ROW_LABELS = ['h277 analogue', r'sech$^2$ sampling']
RFINAL_BINS = [3, 5, 7, 9, 11, 13]
AGE_BINS = [0, 2, 4, 6, 8, 10, 12]
ABSZ_MAX = 3.
DZ = 0.01

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    zfinal_bins = np.arange(-ABSZ_MAX, ABSZ_MAX + DZ, DZ)
    
    # pick discrete colors
    cmap = plt.get_cmap('jet')
    colors = get_color_list(cmap, AGE_BINS)
    
    fig, axs = plt.subplots(len(MIGRATION_SCHEMES), len(RFINAL_BINS)-1, 
                            figsize=(TWO_COLUMN_WIDTH, 3), 
                            sharex=True, sharey=True)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.92, wspace=0., hspace=0.)
    
    for m, row in enumerate(axs):
        mig = MIGRATION_SCHEMES[m]
        output_name = '/'.join((mig, 'insideout/powerlaw_slope11/diskmodel'))
        mzs = MultizoneStars.from_output(output_name)
        row[2].text(0.5, 0.85, ROW_LABELS[m], 
                    ha='center', transform=row[2].transAxes)
        for i, ax in enumerate(row):
            rfinal_lim = tuple(RFINAL_BINS[i:i+2])
            if m == 0:
                ax.set_title(r'$%d\leq R_{\rm{final}} < %d$ kpc' % rfinal_lim)
            for j, color in enumerate(colors):
                age_lim = tuple(AGE_BINS[j:j+2])
                # bin h277 data by formation radius and age
                subset = mzs.filter({'age': age_lim, 'galr_final': rfinal_lim})
                hist, _ = np.histogram(subset('zfinal'), zfinal_bins, 
                                       density=True)
                # apply boxcar smoothing with a width of 0.1 kpc
                hist_smooth = box_smooth(hist, zfinal_bins, 0.1)
                bin_centers = (zfinal_bins[:-1] + zfinal_bins[1:]) / 2
                ax.plot(bin_centers, hist_smooth, c=color, ls='-', 
                        label=r'$%d - %d$ Gyr' % age_lim, zorder=10-j)
    
    axs[0,0].xaxis.set_major_locator(MultipleLocator(1))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_major_locator(MultipleLocator(1))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.2))
    axs[0,0].set_xlim((-2.2, 2.2))
    axs[0,0].set_ylim((0, 4.6))
    for ax in axs[-1]:
        ax.set_xlabel(r'$z_{\rm{final}}$ [kpc]')
    for ax in axs[:,0]:
        ax.set_ylabel('PDF')
    axs[1,-1].legend(loc='upper right', frameon=False, 
                     handlelength=1, handletextpad=0.5)
    plt.savefig(paths.figures / 'midplane_distance.pdf', dpi=300)

if __name__ == '__main__':
    main()
