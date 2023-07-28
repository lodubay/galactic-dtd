"""
This script plots the [O/Fe] distribution within a narrow slice of [Fe/H]
for a VICE multi-zone output, and runs a peak-finding algorithm to determine
whether the output exhibits alpha-element bimodality.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from multizone_stars import MultizoneStars
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

def main(output_name, smoothing=SMOOTH_WIDTH, uncertainties=True):
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig, ax = plt.subplots(figsize=(ONE_COLUMN_WIDTH, ONE_COLUMN_WIDTH), tight_layout=True)
    
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
    
    ax.set_xlabel('[O/Fe]')
    ax.set_ylabel('Normalized PDF')
    ax.set_title(output_name)
    ax.set_xlim((-0.25, 0.55))
    ax.set_ylim((0, None))
    ax.legend(loc='upper left', frameon=False, title='[Fe/H] bin')
    # Save
    fname = output_name.replace('diskmodel', 'ofe_bimodality.png')
    fullpath = paths.figures / 'supplementary' / fname
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    plt.savefig(fullpath, dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_bimodality.py',
        description='Plot the [O/Fe] bimodality test from a multizone run.'
    )
    parser.add_argument('output_name', metavar='NAME',
                        help='Name of VICE multizone output')
    parser.add_argument('--smoothing', metavar='WIDTH', type=float,
                        default=SMOOTH_WIDTH,
                        help='Width of boxcar smoothing (default: 0.05)')
    parser.add_argument('-u', '--uncertainties', action='store_true',
                        help='Model APOGEE uncertainties in VICE output')
    args = parser.parse_args()
    main(**vars(args))
    