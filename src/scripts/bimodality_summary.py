"""
This script plots the [O/Fe] distribution within a narrow slice of [Fe/H]
for a VICE multi-zone output, and runs a peak-finding algorithm to determine
whether the output exhibits alpha-element bimodality.
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from multizone_stars import MultizoneStars
import paths
from utils import get_bin_centers
from _globals import TWO_COLUMN_WIDTH

OFE_LIM = (-0.15, 0.55)
PROMINENCE = 0.2
FEH_BINS = [(-0.6, -0.4), (-0.4, -0.2)]
LINESTYLES = ['--', '-']
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)
SMOOTH_WIDTH = 0.05

SFH_LIST = [
    'insideout', 
    'lateburst', 
    'earlyburst', 
    'twoinfall'
]
DTD_LIST = [
    'powerlaw_slope11', 
    'powerlaw_slope14', 
    'exponential_timescale15',
    'exponential_timescale30', 
    'plateau_width03', 
    'plateau_width10', 
    'prompt',
    'triple'
]

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    figwidth = TWO_COLUMN_WIDTH
    fig, axs = plt.subplots(len(DTD_LIST), len(SFH_LIST), 
                            figsize=(figwidth, figwidth*1.67), 
                            sharex=True, sharey=True)
    fig.subplots_adjust(top=0.98, bottom=0.03, left=0.07, right=0.98, wspace=0., hspace=0.)
    
    with tqdm(total=len(SFH_LIST) * len(DTD_LIST)) as t:
        for j, evolution in enumerate(SFH_LIST):
            axs[0,j].set_title(evolution)
            for i, RIa in enumerate(DTD_LIST):
                axs[i,0].set_ylabel(RIa)
                # Import VICE multi-zone output data
                output_name = '/'.join(['gaussian', evolution, RIa, 'diskmodel'])
                plot_bimodality(axs[i,j], output_name, uncertainties=True)
                t.update()
    
    for ax in axs[-1]:
        ax.set_xlabel('[O/Fe]')
    axs[0,0].set_xlim((-0.25, 0.55))
    axs[0,0].set_ylim((0, 1.1))
    axs[0,0].legend(loc='upper left', frameon=False, title='[Fe/H] bin')
    # Save
    plt.savefig(paths.figures / 'bimodality_summary.png', dpi=300)
    plt.close()
    

def plot_bimodality(ax, output_name, smoothing=SMOOTH_WIDTH, uncertainties=True):
    mzs = MultizoneStars.from_output(output_name)
    if uncertainties:
        mzs.model_uncertainty(inplace=True)
    subset = mzs.region(galr_lim=GALR_LIM, absz_lim=ABSZ_LIM)
    for i, feh_bin in enumerate(FEH_BINS):
        subset_slice = subset.filter({'[fe/h]': feh_bin})
        mdf, bin_edges = subset_slice.mdf('[o/fe]', smoothing=smoothing,
                                          bins=np.arange(-0.15, 0.56, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        peaks, _ = find_peaks(mdf, prominence=PROMINENCE)
        ax.plot(bin_centers, mdf, ls=LINESTYLES[i], label=feh_bin)
        ax.plot(bin_centers[peaks], mdf[peaks], 'rx')


if __name__ == '__main__':
    main()
    