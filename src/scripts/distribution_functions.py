"""
This file contains generic functions for plotting distributions of data
binned by galactic radius and z-height.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
from utils import get_bin_centers, discrete_colormap, get_color_list
from _globals import ABSZ_BINS, GALR_BINS

# Ratio between major and minor tick spacing in data units
MAJOR_MINOR_RATIO = 5.

    
def plot_distributions(func, data, axs, label='', cmap_name='plasma_r',
                       galr_bins=GALR_BINS, absz_bins=ABSZ_BINS):
    """
    Plot distributions of data binned by galactic radius and z-height.
    
    Parameters
    ----------
    func : <function>
        The function which computes the distribution in the given bin of
        galactic coordinates. Must take data as an argument and bounds on
        galr and absz as keyword arguments, and must return the values of the
        distribution and its bin edges (similar to the output of 
        numpy.histogram).
    data : pandas.DataFrame
        Data (e.g. from VICE or APOGEE) from which to generate the 
        distributions.
    axs : list of matplotlib.axes.Axes
        Axes on which to plot the age distributions, the length of which must
        correspond to len(absz_bins)-1 (3 by default); usually a single
        column from a larger array of axes.
    label : str, optional
        Axis column label. The default is 'VICE'.
    cmap_name : str, optional
        The name of the colormap to code by galactic radius. The default is
        'plasma_r'.
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
    galr_bins : list, optional
        Bin edges of galactic radius in kpc. The default is 
        [3, 5, 7, 9, 11, 13, 15].
    """
    cmap = plt.get_cmap(cmap_name)
    colors = get_color_list(cmap, galr_bins)
    if len(colors) != len(galr_bins) - 1:
        colors = [None for galr in galr_bins[:-1]]
    if len(axs) == len(absz_bins) - 1:
        for i, ax in enumerate(axs.flatten()):
            absz_lim = absz_bins[-(i+2):len(absz_bins)-i]
            for j in range(len(galr_bins)-1):
                galr_lim = galr_bins[j:j+2]
                dist, bin_edges = func(data, galr_lim, absz_lim)
                ax.plot(get_bin_centers(bin_edges), dist, 
                        color=colors[j], linewidth=1)
        axs[0].set_title(label)
    else:
        raise ValueError('Mismatch between axes and z-height bins.')
    

def setup_axes(ncols=2, figure_width=3.25, xlabel='', xlim=None, 
               major_tick_spacing=1, galr_bins=GALR_BINS, absz_bins=ABSZ_BINS, 
               cbar_width=0.6, cmap_name='plasma_r', panel_aspect_ratio=1.5):
    """
    Set up matplotlib figure and axes for the distribution plot.
    
    Parameters
    ----------
    ncols : int, optional
        Number of columns in figure. The default is 2.
    figure_width : float, optional
        Width of the figure in inches. For an AASTeX document, this should be
        3.25 for a single-column figure or 7.0 for a double-column figure. The
        figure height will be scaled to maintain a consistent aspect ratio.
        The default is 3.25.
    xlabel : str, optional
        The x-axis label. The default is ''.
    xlim : tuple or NoneType, optional
        The x-axis limits. The default is None.
    major_tick_spacing : float, optional
        The spacing in data units between major tick marks. The default is 1.
    galr_bins : list, optional
        Bin edges of galactic radius in kpc. The default is 
        [3, 5, 7, 9, 11, 13, 15].
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
    cbar_width : float, optional
        The width of the colorbar as a fraction of the figure size. The default
        is 0.6.
    panel_aspect_ratio : float, optional
        The aspect ratio of the panels, equal to the panel width divided by
        the panel height. The default is 1.5.
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : list of matplotlib.axes.Axes
    """
    # Determine figure dimensions
    nrows = len(absz_bins) - 1
    ax_width = (figure_width - 0.25) / ncols
    ax_height = ax_width / panel_aspect_ratio
    figure_height = ax_height * nrows + 1
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey='row',
                            figsize=(figure_width, figure_height))
    fig.subplots_adjust(left=0.12, top=0.93, right=0.97, bottom=0.22,
                        wspace=0.07, hspace=0.07)
    # Format x-axis
    axs[0,0].set_xlim(xlim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(major_tick_spacing))
    minor_tick_spacing = major_tick_spacing / MAJOR_MINOR_RATIO
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(minor_tick_spacing))
    for ax in axs[-1]:
        ax.set_xlabel(xlabel)
    # Remove spines and y-axis labels
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
        # Set bottom ticks pointing out
        ax.tick_params(axis='x', which='both', direction='out')
    # Add common y-axis label
    fig.text(0.03, 0.58, r'Distance from Galactic midplane $|z|$',
              ha='center', va='center', rotation='vertical')
    # Label rows
    for i in range(len(absz_bins)-1):
        absz_lim = tuple(absz_bins[-(i+2):len(absz_bins)-i])
        axs[i,0].set_ylabel(r'$%s - %s$ kpc' % absz_lim)
    # Add colorbar on bottom
    cmap, norm = discrete_colormap(cmap_name, galr_bins)
    cax = plt.axes([0.5 - (cbar_width / 2), 0.09, cbar_width, 0.02])
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax,
                        orientation='horizontal')
    cbar.set_label('Galactocentric radius [kpc]')
    return fig, axs
