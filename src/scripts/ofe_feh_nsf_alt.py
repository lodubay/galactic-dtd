"""
This script plots [O/Fe] vs [Fe/H] for several VICE multizone simulation outputs.
Intended for use in an NSF grant proposal.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MultipleLocator
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    sample_dataframe
from _globals import END_TIME, ZONE_WIDTH

OUTPUTS = ["diffusion/conroy22/powerlaw_slope11",
           "diffusion/spitoni21/powerlaw_slope11",
           "post-process/spitoni21/powerlaw_slope11"]
DATA_DIR = paths.data / "migration"
GALR_BINS = [(4, 6), (8, 10), (12, 15)]
FIGSIZE = (5, 4.5)
FEH_LIM = (-2.6, 1)
OFE_LIM_LIST = [(-0.1, 0.75), (-0.2, 0.5), (-0.2, 0.5)]
CMAP = "plasma_r"
LABELS = ["Exponential infall w/ Conroy+ 2022 SFE", "Two-infall w/ mixing", "Two-infall no mixing"]
GALR_MAX = 15.5


def main():
    fig, axs = setup_axes(figsize=FIGSIZE, xlim=FEH_LIM, ylim_list=OFE_LIM_LIST)
    norm = Normalize(vmin=-2, vmax=GALR_MAX)
    cax = setup_colorbar(fig, CMAP, norm, label="Birth radius [kpc]")
    
    for i, output in enumerate(OUTPUTS):
        stars = multioutput_to_pandas(output, DATA_DIR)
        plot_stars(stars, axs[i,:], galr_bins=GALR_BINS, cmap=CMAP, norm=norm,
                   zone_origin=("post-process" in output))
                   
    # SFH labels
    for ax, label in zip(axs[:,1], LABELS):
        ax.set_title(label)
    plt.savefig(paths.figures / "ofe_feh_nsf_alt.png", dpi=300)
    plt.close()
        
        
def plot_stars(stars, axs, galr_bins=[(4, 6), (8, 10), (12, 15)], cmap="winter",
               norm=None, zone_origin=False):
    """
    Plot [O/Fe] vs [Fe/H] color-coded by age for the given star particles.
    
    Parameters
    ----------
    stars : pandas.DataFrame
        Star particle data from a VICE multizone output.
    axs : list of matplotlib.axes.Axes
        Single row of axes on which to plot the stars data.
    galr_bins : list of tuples
        Bins of galactocentric radius in which to include stars data.
    """
    for (ax, galr_lim) in zip(axs, galr_bins):
        filtered = filter_multioutput_stars(stars, galr_lim=galr_lim, 
                                            zone_origin=zone_origin)
        sample = sample_dataframe(filtered, 10000)
        ax.scatter(sample["[fe/h]"], sample["[o/fe]"], s=0.5,
                   c=sample["zone_origin"] * ZONE_WIDTH, cmap=cmap, norm=norm,
                   rasterized=True, edgecolor="none")


def setup_colorbar(fig, cmap, norm, label=""):
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
    plt.subplots_adjust(right=0.85, left=0.1, bottom=0.07, top=0.95,
                        wspace=0.05, hspace=0.25)
    cax = plt.axes([0.88, 0.17, 0.03, 0.68])
    # Add colorbar
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax)
    cbar.set_label(label)
    cax.yaxis.set_minor_locator(MultipleLocator(0.5))
    cax.set_ylim((0, None))
    return cax


def setup_axes(figsize=(4, 4), xlim=(-3, 1), 
               ylim_list=[(-0.1, 0.6), (-0.1, 0.65), (-0.2, 0.5), (-0.2, 0.5)],
               galr_bins=[(4, 6), (8, 10), (12, 15)]):
    """
    Set up figure and grid of axes.
    
    Parameters
    ----------
    figsize : tuple, optional
        Figure size (width, height) in inches. The default is (4, 4).
    """
    nrows = len(OUTPUTS)
    ncols = len(galr_bins)
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey="row")
    # Axis labels and limits
    for ax in axs[-1,:]:
        ax.set_xlabel("[Fe/H]")
        ax.set_xlim(xlim)
    for ax, ylim in zip(axs[:,0], ylim_list): 
        ax.set_ylabel("[O/Fe]")
        ax.set_ylim(ylim)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    # Axis titles
    for (ax, galr_lim) in zip(axs[0,:], galr_bins):
        galr_label = r"%s kpc $\leq R_{\rm gal} <$ %s kpc" % galr_lim
        ax.text(-2.4, 0.7, galr_label, size=8, va="top")
    # Tick marks
    axs[0,0].xaxis.set_major_locator(MultipleLocator(1.))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.2))
    return fig, axs


if __name__ == "__main__":
    main()
