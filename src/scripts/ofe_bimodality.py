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
from apogee_tools import import_apogee, apogee_region, apogee_mdf
from colormaps import paultol
from _globals import TWO_COLUMN_WIDTH, ONE_COLUMN_WIDTH

OFE_LIM = (-0.15, 0.55)
PROMINENCE = 0.1
FEH_BINS = [(-0.6, -0.4), (-0.4, -0.2)]
LINESTYLES = ['-', '--']
COLORS = [paultol.highcontrast.colors[2], paultol.highcontrast.colors[0]]
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)
SMOOTH_WIDTH = 0.05

def main(output_name, smoothing=SMOOTH_WIDTH, uncertainties=True, **kwargs):
    apogee_data = import_apogee()
    mzs = MultizoneStars.from_output(output_name)
    if uncertainties:
        mzs.model_uncertainty(inplace=True, apogee_data=apogee_data)
    plot_bimodality_comparison(mzs, apogee_data, **kwargs)
    
    
def plot_bimodality_comparison(mzs, apogee_data, style='paper', **kwargs):
    """
    Plot the [O/Fe] bimodality from a VICE output and from APOGEE.
    
    Parameters
    ----------
    mzs : MultizoneStars object
    apogee_data : pandas DataFrame
    style : str, optional
        Plotting style to use. The default is 'paper'.
    **kwargs passed to plot_bimodality()
    """
    plt.style.use(paths.styles / f'{style}.mplstyle')
    fig, axs = plt.subplots(1, 2, figsize=(TWO_COLUMN_WIDTH, ONE_COLUMN_WIDTH), 
                            tight_layout=True, sharex=True, sharey=True)
    # Plot VICE bimodality
    plot_bimodality(axs[0], mzs, apogee_data=apogee_data, **kwargs)
    # Plot APOGEE bimodality
    subset = apogee_region(apogee_data, galr_lim=GALR_LIM, absz_lim=ABSZ_LIM)
    for i, feh_bin in enumerate(FEH_BINS):
        subset_slice = subset[(subset['FE_H'] >= feh_bin[0]) & (subset['FE_H'] < feh_bin[1])]
        mdf, bin_edges = apogee_mdf(subset_slice, col='O_FE', 
                                    smoothing=SMOOTH_WIDTH, 
                                    bins=np.arange(OFE_LIM[0], OFE_LIM[1]+0.01, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        peaks, _ = find_peaks(mdf, prominence=PROMINENCE)
        axs[1].plot(bin_centers, mdf, ls=LINESTYLES[i], c=COLORS[i], 
                    label=feh_bin)
        axs[1].plot(bin_centers[peaks], mdf[peaks], 'rx')
    
    axs[0].set_xlabel('[O/Fe]')
    axs[1].set_xlabel('[O/Fe]')
    axs[0].set_ylabel('Normalized PDF')
    fig.suptitle(mzs.name)
    axs[0].set_xlim(OFE_LIM)
    axs[0].set_ylim((0, None))
    axs[0].legend(loc='best', frameon=False, title='[Fe/H] bin')
    # Save
    fname = mzs.name.replace('diskmodel', 'ofe_bimodality.png')
    fullpath = paths.extra / fname
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    plt.savefig(fullpath, dpi=300)
    plt.close()


def plot_bimodality(ax, mzs, feh_bins=FEH_BINS, smoothing=SMOOTH_WIDTH, 
                    resample=True,  apogee_data=None, 
                    nsamples=100000, ofe_lim=OFE_LIM, galr_lim=GALR_LIM,
                    absz_lim=ABSZ_LIM, linestyles=LINESTYLES, colors=COLORS,
                    show_peaks=True, prominence=0.1, **kwargs):
    """
    Plot the [O/Fe] distribution for the given slices of [Fe/H].
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object on which to plot the distributions.
    mzs : MultizoneStars object
        Instance with VICE multi-zone star output
    feh_bins : list of tuples, optional
        Bins in [Fe/H] which define the slices for the [O/Fe] distributions.
        The default is [(-0.6, -0.4), (-0.4, -0.2)].
    smoothing : float, optional
        Smoothing width for the [O/Fe] distribution. The default is 0.05.
    resample : bool, optional
        If True, resample the |z| distribution of stars to match APOGEE. If
        False, all model stellar populations are counted. The default is True.
    apogee_data : pandas.DataFrame or NoneType, optional
        DataFrame containing APOGEE abundance data. If resample==True or 
        uncertainties==True and apogee_data==None, the APOGEE data file will
        be imported. The default is None.
    nsamples : int, optional
        Number of stellar populations to sample if resample==True. The default
        is 20000.
    ofe_lim : tuple, optional
        Minimum and maximum [O/Fe] to be included in the distribution.
        The default is (-0.15, 0.55)
    galr_lim : tuple, optional
        Limits on galactocentric radius. The default is (7, 9) kpc.
    absz_lim : tuple, optional
        Limits on absolute midplane distance. The default is (0, 2) kpc.
    colors : list, optional
        List of line colors.
    linestyles : list, optional
        List of line styles.
    show_peaks : bool, optional
        Whether to indicate the location of prominent peaks. The default is
        True.
    prominence : float, optional
        Prominence threshold for what counts as a peak. The default is 0.1.
        
    **kwargs passed to ax.plot
    
    Returns
    -------
    list of matplotlib.line.Line2D
    """
    # Import APOGEE data if needed
    if resample and apogee_data is None:
        apogee_data = import_apogee()
    subset = mzs.region(galr_lim=galr_lim, absz_lim=absz_lim)
    if resample:
        apogee_subset = apogee_region(apogee_data, galr_lim=galr_lim, 
                                      absz_lim=absz_lim)
        subset.resample_zheight(nsamples, apogee_subset, inplace=True)
    lines = []
    for i, feh_bin in enumerate(feh_bins):
        subset_slice = subset.filter({'[fe/h]': feh_bin})
        mdf, bin_edges = subset_slice.mdf('[o/fe]', smoothing=smoothing,
                                          bins=np.arange(ofe_lim[0], 
                                                         ofe_lim[1]+0.01, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        line = ax.plot(bin_centers, mdf, label=feh_bin, c=colors[i], 
                       ls=linestyles[i])
        lines.append(line)
        if show_peaks:
            peaks, _ = find_peaks(mdf, prominence=prominence)
            ax.plot(bin_centers[peaks], mdf[peaks], 'rx')
    return lines


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
    parser.add_argument('-r', '--resample', action='store_true',
                        help='Re-sample |z| distribution of model to match APOGEE')
    args = parser.parse_args()
    main(**vars(args))
    