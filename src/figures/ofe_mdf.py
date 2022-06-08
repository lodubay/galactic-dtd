"""
Plot metallicity distribution functions (MDFs) of [O/Fe] binned by radius.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vice
from ofe_feh_vice import GALR_BINS, ABSZ_BINS, FEH_LIM, OFE_LIM, ZONE_WIDTH
from ofe_feh_apogee import apogee_region
from utils import multioutput_to_pandas, filter_multioutput_stars

def main(output_name, migration_dir='../data/migration_outputs',
         apogee_path='../data/APOGEE/dr17_cut_data.csv', cmap_name='copper'):
    print('Importing data')
    stars = multioutput_to_pandas(output_name, migration_dir)
    data = pd.read_csv(Path(apogee_path))

    fig, axs = plt.subplots(3, 2, figsize=(6, 9), sharex=True)
    cmap = plt.get_cmap(cmap_name)
    colors = cmap([r/15 for r in GALR_BINS[:-1]])
    for i in range(len(ABSZ_BINS)-1):
        print('Panel %s' % (i+1))
        absz_lim = ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i]
        for j in range(len(GALR_BINS)-1):
            galr_lim = GALR_BINS[j:j+2]
            # Plot VICE in left panels
            # vice_subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
            #                                        ZONE_WIDTH)
            # axs[i,0].hist(vice_subset['[o/fe]'], bins=80, range=(-0.2, 0.6),
            #             density=True, histtype='step', color=colors[j])
            mdf_sum = 0
            for z in range(int(galr_lim[0] / ZONE_WIDTH), int(galr_lim[1] / ZONE_WIDTH)):
                singlezone_output = Path('%s.vice/zone%s' % (output_name, z))
                mdf = vice.mdf(str(Path(migration_dir) / singlezone_output))
                mdf_sum += np.array(mdf['dn/d[o/fe]'])
            bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
            axs[i,0].hist(bins[:-1], bins, weights=mdf_sum, histtype='step',
                          color=colors[j], density=True)

            # Plot APOGEE in right panels
            apogee_subset = apogee_region(data, galr_lim, absz_lim)
            axs[i,1].hist(apogee_subset['O_FE'], bins=80, range=(-0.2, 0.6),
                        density=True, histtype='step', color=colors[j])
    axs[0,0].set_xlim((-0.2, 0.6))
    plt.savefig('ofe_mdf.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main('diffusion/insideout/powerlaw')