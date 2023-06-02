"""
This script plots the distributions of migration distance as a function of
stellar population age and formation radius.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from utils import multioutput_to_pandas, box_smooth, get_color_list
import paths

evolution = 'insideout'
RIa = 'powerlaw_slope11'
migration_schemes = ['diffusion', 'gaussian']
labels = ['h277 analog', 'Gaussian sampling']

rform_bins = np.arange(3, 15, 2, dtype='float')
age_bins = np.arange(0, 14, 2, dtype='float')
zfinal_bins = np.arange(-5, 5.01, 0.01)

# pick discrete colors
cmap = plt.get_cmap('jet')
colors = get_color_list(cmap, age_bins)

fig, axs = plt.subplots(2, 5, figsize=(15, 6), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0., hspace=0.)

for m, row in enumerate(axs):
    mig = migration_schemes[m]
    output_name = '/'.join((mig, evolution, RIa))
    stars = multioutput_to_pandas(output_name)
    row[2].text(0.9, 0.9, labels[m], ha='right', transform=row[2].transAxes)
    for i, ax in enumerate(row):
        rform_lim = tuple(rform_bins[i:i+2])
        ax.vlines(rform_lim, ymin=0, ymax=0.6, linestyle=':', color='k')
        if m == 0:
            ax.set_title(r'$%d\leq R_{\rm{form}} < %d$ kpc' % rform_lim)
        for j, color in enumerate(colors):
            age_lim = tuple(age_bins[j:j+2])
            # bin h277 data by formation radius and age
            subset = stars[(stars['galr_origin'] >= rform_lim[0]) &
                           (stars['galr_origin'] <  rform_lim[1]) &
                           (stars['age']         >= age_lim[0]) &
                           (stars['age']         <  age_lim[1])]
            # limit to bins with a meaningful number of stars
            no_dups = subset.drop_duplicates('zfinal')
            if no_dups.shape[0] > 250:
                hist, _ = np.histogram(subset['zfinal'], zfinal_bins, 
                                       density=True)
                # apply boxcar smoothing with a width of 0.1 kpc
                hist_smooth = box_smooth(hist, zfinal_bins, 0.1)
                bin_centers = (zfinal_bins[:-1] + zfinal_bins[1:]) / 2
                ax.plot(bin_centers, hist_smooth, c=color, ls='-', 
                        label=r'$%d - %d$ Gyr' % age_lim)

axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
axs[0,0].set_xlim((-2.1, 2.1))
axs[0,0].set_ylim((0, 4.0))
for ax in axs[-1]:
    # ax.set_xlabel(r'$\Delta R_{\rm{gal}}$ [kpc]')
    ax.set_xlabel(r'$z_{\rm{final}}$ [kpc]')
for ax in axs[:,0]:
    ax.set_ylabel('PDF')
axs[1,-1].legend(title='Age bins', loc='upper right', fontsize=8, frameon=False)
plt.savefig(paths.figures / 'check_vertical_migration.png', dpi=300)
