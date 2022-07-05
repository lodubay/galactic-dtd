"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights.
"""

import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import vice
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    sample_dataframe

global FEH_LIM
global OFE_LIM
global GALR_BINS
global ABSZ_BINS
global ZONE_WIDTH

FEH_LIM = (-1.3, 0.8)
OFE_LIM = (-0.1, 0.5)
GALR_BINS = [3, 5, 7, 9, 11, 13] # kpc
ABSZ_BINS = [0, 0.5, 1, 2] # kpc
ZONE_WIDTH = 0.1 # kpc

def main(output_name, data_dir='../data/migration_outputs', cmap='winter'):
    """
    Parameters
    ----------
    output_name : string
        Path to the
    """
    # Import multioutput stars data
    stars = multioutput_to_pandas(output_name, data_dir)
    fig, axs = plot_ofe_feh_stars(stars, cmap)
    plot_post_process_track(output_name, axs, galr=8, data_dir=data_dir)
    plot_post_process_tracks(output_name, axs, data_dir=data_dir)
    plt.savefig('ofe_feh_vice.pdf', dpi=300)
    plt.close()


def plot_ofe_feh_stars(stars, cmap):
    """
    Plot just the VICE multizone stars data
    """
    fig, axs = setup_axes(len(ABSZ_BINS)-1, len(GALR_BINS)-1,
                               xlim=FEH_LIM, ylim=OFE_LIM)
    norm = normalize_colorbar(stars)
    cax = setup_colorbar(fig, cmap, norm, label=r'Birth $R_{\rm{Gal}}$ [kpc]')
    cax.yaxis.set_minor_locator(MultipleLocator(0.5))
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = filter_multioutput_stars(stars, galr_lim, absz_lim,
                                              ZONE_WIDTH)
            sample = sample_dataframe(subset, 10000)
            # Scatter plot of random sample of stellar particles
            ax.scatter(sample['[fe/h]'], sample['[o/fe]'], s=0.1,
                       c=sample['zone_origin'] * ZONE_WIDTH, cmap=cmap,
                       norm=norm, rasterized=True)
            # Label axes
            if i == len(axs)-1:
                ax.set_xlabel('[Fe/H]')
            if j == 0:
                ax.set_ylabel('[O/Fe]')
                ax.text(0.55, 0.85, r'$%s\leq |z| < %s$' % absz_lim,
                        transform=ax.transAxes, size=8)
            if i == 0:
                ax.set_title(r'$%s\leq R_{\rm{Gal}} < %s$ kpc'% galr_lim)
    return fig, axs


def plot_post_process_track(output_name, axs, galr=8,
                            data_dir='../data/migration_outputs'):
    """
    Plot abundance track for a given annulus across all panels.

    Parameters
    ----------
    output_name : str
        Path to the VICE output directory
    axs : list of axes
    galr : float, optional [default: 8]
        Galactic radius of abundance track to plot in kpc. The default is the
        solar annulus
    data_dir : str, optional
        Parent directory for all migration outputs
    """
    # Convert annulus from kpc to zone number
    zone = int(galr / ZONE_WIDTH)
    # Import post-processed output for the given annulus
    post_process_path = Path(data_dir) / Path(output_name.replace(
        'diffusion', 'post-process') + '.vice') / Path('zone%d' % zone)
    post_process_hist = vice.history(str(post_process_path))
    # Plot abundance track on all panels
    for ax in axs.flatten():
        ax.plot(post_process_hist['[fe/h]'], post_process_hist['[o/fe]'],
                c='k', ls='--', linewidth=1)


def plot_post_process_tracks(output_name, axs,
                             data_dir='../data/migration_outputs'):
    """
    Plot post-process abundance track for the mean annulus in each panel.

    Parameters
    ----------
    output_name : str
        Path to the VICE output directory
    axs : list of axes
    data_dir : str, optional
        Parent directory for all migration outputs
    """
    # Convert annulus from kpc to zone number
    nzones = len(GALR_BINS)-1
    mean_galr = [(GALR_BINS[i] + GALR_BINS[i+1]) / 2 for i in range(nzones)]
    zones = [int(mean_galr[i] / ZONE_WIDTH) for i in range(nzones)]
    for i, zone in enumerate(zones):
        # Import post-processed output for the given annulus
        post_process_path = Path(data_dir) / Path(output_name.replace(
            'diffusion', 'post-process') + '.vice') / Path('zone%d' % zone)
        post_process_hist = vice.history(str(post_process_path))
        for ax in axs[:,i].flatten():
            ax.plot(post_process_hist['[fe/h]'], post_process_hist['[o/fe]'],
                    c='k', ls='-', linewidth=1)


def normalize_colorbar(stars):
    """
    Calculate the colorbar normalization for the full migration output data
    before subsampling.

    Parameters
    ----------
    stars : pandas DataFrame
        VICE migration output

    Returns
    -------
    norm : instance of matplotlib.colors.Normalize
    """
    # Filter to full range to be plotted
    subset = filter_multioutput_stars(stars,
                                      galr_lim=(GALR_BINS[0], GALR_BINS[-1]),
                                      absz_lim=(ABSZ_BINS[0], ABSZ_BINS[-1]),
                                      zone_width=ZONE_WIDTH)
    norm = Normalize(vmin=0, vmax=subset['zone_origin'].max() * ZONE_WIDTH)
    return norm


def setup_colorbar(fig, cmap, norm, label=''):
    """
    Configure colorbar given a colormap and a normalization.

    Parameters
    ----------
    fig
    cmap
    norm

    Returns
    -------
    cax : colorbar axis
    """
    # Colorbar axis
    plt.subplots_adjust(right=0.92, left=0.06, bottom=0.07, top=0.95,
                        wspace=0.05, hspace=0.05)
    cax = plt.axes([0.93, 0.07, 0.02, 0.88])
    # Add colorbar
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label(label)
    return cax


def setup_axes(rows, cols, width=8, xlim=None, ylim=None):
    """
    Set up a blank grid of axes plus a colorbar axis.

    Parameters
    ----------
    rows : int
        Number of rows of axes
    cols : int
        Number of columns of axes
    width : float, optional
        Width of the figure in inches. The default is 8 in.
    xlim : tuple or None, optional
        Limits of x-axis for all axes
    ylim : tuple or None, optional
        Limits of y-axis for all axes

    Returns
    -------
    fig : matplotlib figure
    axs : list of axes
    cax : axis object for colorbar
    """
    fig, axs = plt.subplots(rows, cols, figsize=(width, (width/cols)*rows),
                            sharex=True, sharey=True)
    # Configure axis limits and ticks (will be applied to all axes)
    axs[0,0].set_xlim(xlim)
    axs[0,0].set_ylim(ylim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.02))
    return fig, axs


if __name__ == '__main__':
    # temporary sys arguments - replace with argparse later
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))