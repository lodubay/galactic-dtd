"""
This script plots the [O/Fe] distribution within a narrow slice of [Fe/H]
for APOGEE data, and runs a peak-finding algorithm to determine
whether the output exhibits alpha-element bimodality.
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from apogee_tools import import_apogee, apogee_region, apogee_mdf
import paths
from utils import get_bin_centers
from _globals import ONE_COLUMN_WIDTH

OFE_LIM = (-0.15, 0.55)
PROMINENCE = 0.1
FEH_BINS = [(-0.6, -0.4), (-0.4, -0.2)]
LINESTYLES = ['--', '-']
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)
SMOOTH_WIDTH = 0.05

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    figwidth = ONE_COLUMN_WIDTH
    fig, ax = plt.subplots(figsize=(figwidth, figwidth), 
                           tight_layout=True)
    
    data = import_apogee()
    subset = apogee_region(data, galr_lim=GALR_LIM, absz_lim=ABSZ_LIM)
    for i, feh_bin in enumerate(FEH_BINS):
        subset_slice = subset[(subset['FE_H'] >= feh_bin[0]) & (subset['FE_H'] < feh_bin[1])]
        mdf, bin_edges = apogee_mdf(subset_slice, col='O_FE', 
                                    smoothing=SMOOTH_WIDTH, 
                                    bins=np.arange(-0.15, 0.56, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        peaks, _ = find_peaks(mdf, prominence=PROMINENCE)
        ax.plot(bin_centers, mdf, ls=LINESTYLES[i], label=feh_bin)
        ax.plot(bin_centers[peaks], mdf[peaks], 'rx')
    
    ax.set_xlabel('[O/Fe]')
    ax.set_ylabel('Normalized PDF')
    ax.set_title('APOGEE DR17')
    ax.set_xlim((-0.25, 0.55))
    ax.set_ylim((0, 1.1))
    ax.legend(loc='upper left', frameon=False, title='[Fe/H] bin')
    # Save
    plt.savefig(paths.figures / 'bimodality_apogee.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
    