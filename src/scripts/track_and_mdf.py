"""
This file contains functions for plotting abundance tracks in [Fe/H] vs [O/Fe]
alongside their corresponding metallicity distribution functions (MDFs).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# from utils import get_bin_centers


def plot_vice_onezone(output, fig=None, axs=[], label=None, color=None,
                      histtype='step', logmdf=False, marker_labels=False,
                      style_kw={}):
    """
    Wrapper for plot_track_and_mdf given a VICE onezone output.

    Parameters
    ----------
    output : str
        Path to VICE output, not necessarily including the '.vice' extension
    fig : instance of matplotlib.figure.figure, optional
        If no figure is provided, one is generated from setup_axes
    axs : list of matplotlib.axes.Axes
        There should be three axes: the first for the main [Fe/H] vs [O/Fe]
        panel, the second for the MDF in [Fe/H], and the third for the MDF
        in [O/Fe]
    label : str, optional
        Plot label to add to main panel legend
    color : str, optional
        Color of plot and histogram
    histtype : str, optional
        Histogram type; options are 'bar', 'barstacked', 'step', 'stepfilled'.
        The default is 'step'.
    logmdf : bool, optional
        If True, plot the marginal distributions on a log scale. The default
        is False.
    style_kw : dict, optional
        Dict of style-related keyword arguments to pass to both
        matplotlib.pyplot.plot and matplotlib.pyplot.hist

    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : list of matplotlib.axes.Axes
    """
    import vice
    hist = vice.history(output)
    mdf = vice.mdf(output)
    mdf_bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
    fig, axs = plot_track_and_mdf(hist['[fe/h]'], hist['[o/fe]'],
                                  dn_dfeh=mdf['dn/d[fe/h]'], feh_bins=mdf_bins,
                                  dn_dofe=mdf['dn/d[o/fe]'], ofe_bins=mdf_bins,
                                  fig=fig, axs=axs, label=label, color=color,
                                  histtype=histtype, logmdf=logmdf,
                                  style_kw=style_kw)
    if color == None:
        color = axs[0].lines[-1].get_color()
    plot_time_markers(hist['time'], hist['[fe/h]'], hist['[o/fe]'], axs[0],
                      color=color, show_labels=marker_labels)
    return fig, axs


def plot_time_markers(time, feh, ofe, ax, loc=[0.1, 0.3, 1, 3, 10],
                      color=None, show_labels=False):
    """
    Add temporal markers to the [O/Fe] vs [Fe/H] tracks.

    Parameters
    ----------
    time : array-like
        Array of simulation time in Gyr.
    feh : array-like
        Array of [Fe/H] abundances.
    ofe : array-like
        Array of [O/Fe] abundances.
    ax : matplotlib.axes.Axes
        Axis in which to plot the time markers.
    loc : list, optional
        List of times in Gyr to mark. The default is [0.1, 0.3, 1, 3, 10].
    color : color or None, optional
        Color of the markers. The default is None.
    show_labels : bool, optional
        Whether to add marker labels. The default is False.
    """
    markers = ['o', 's', '^', 'd', 'v', 'p', '*', 'X']
    time = np.array(time)
    for i, t in enumerate(loc):
        idx = np.argmin(np.abs(time - t))
        ax.scatter(feh[idx], ofe[idx], s=9, marker=markers[i],
                   edgecolors=color, facecolors='w', zorder=10)
        if show_labels:
            if t < 1:
                label = f'{int(t*1000)} Myr'
            else:
                label = f'{int(t)} Gyr'
            if i == 0:
                xpad = -0.1
                ypad = 0.01
            else:
                xpad = 0.03
                ypad = 0.008
            ax.text(feh[idx] + xpad, ofe[idx] + ypad, label, fontsize=7,
                    ha='left', va='bottom', zorder=10,
                    bbox={
                        'facecolor': 'w',
                        'edgecolor': 'none',
                        'boxstyle': 'round',
                        'pad': 0.05,
                        'alpha': 0.8,
                    },
            )


def plot_track_and_mdf(feh, ofe, dn_dfeh=[], feh_bins=10, dn_dofe=[],
                       ofe_bins=10, fig=None, axs=[], label=None, color=None,
                       histtype='step', logmdf=False, style_kw={}):
    """
    Simultaneously plot a track in [Fe/H] vs [O/Fe] and the corresponding
    metallicity distribution functions (MDFs).

    Parameters
    ----------
    feh : array-like
        Array of [Fe/H] abundances
    ofe : array-like
        Array of [O/Fe] abundances
    dn_dfeh : array-like, optional
        Distribution of metallicities dN/d[Fe/H]. If empty, the feh parameter
        is passed to matplotlib.pyplot.hist directly
    feh_bins : int or array-like, optional
        Bins in [Fe/H] for plotting the MDF in [Fe/H]
    dn_dofe : array-like, optional
        Distribution of O abundances dN/d[O/Fe]. If empty, the ofe parameter
        is passed to matplotlib.pyplot.hist directly
    ofe_bins : int or array-like, optional
        Bins in [O/Fe] for plotting the MDF in [O/Fe]
    fig : instance of matplotlib.figure.figure, optional
        If no figure is provided, one is generated from setup_axes
    axs : list of matplotlib.axes.Axes
        There should be three axes: the first for the main [Fe/H] vs [O/Fe]
        panel, the second for the MDF in [Fe/H], and the third for the MDF
        in [O/Fe]
    label : str, optional
        Plot label to add to main panel legend
    color : str, optional
        Color of plot and histogram
    histtype : str, optional
        Histogram type; options are 'bar', 'barstacked', 'step', 'stepfilled'.
        The default is 'step'.
    logmdf : bool, optional
        If True, plot the marginal distributions on a log scale. The default
        is False.
    style_kw : dict, optional
        Dict of style-related keyword arguments to pass to both
        matplotlib.pyplot.plot and matplotlib.pyplot.hist

    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : list of matplotlib.axes.Axes
    """
    if fig == None or len(axs) != 3:
        fig, axs = setup_axes(logmdf=logmdf)

    # Plot abundance tracks on main panel
    axs[0].plot(feh, ofe, label=label, color=color, **style_kw)

    # Plot distribution of [Fe/H] on top panel
    if len(dn_dfeh) == 0:
        dn_dfeh, feh_bins = np.histogram(feh, bins=feh_bins)
    # mask zeros before taking log
    dn_dfeh = np.array(dn_dfeh)
    if logmdf:
        dn_dfeh[dn_dfeh == 0] = 1e-10
        weights = np.log10(dn_dfeh)
    else:
        weights= dn_dfeh
    axs[1].hist(feh_bins[:-1], feh_bins, weights=weights, color=color,
                histtype=histtype, **style_kw)
    # feh_bin_centers = get_bin_centers(feh_bins)
    # axs[1].plot(feh_bin_centers, weights, color=color, **style_kw)

    # Plot distribution of [O/Fe] on right side panel
    if len(dn_dofe) == 0:
        dn_dofe, ofe_bins = np.histogram(ofe, bins=ofe_bins)
    dn_dofe = np.array(dn_dofe)
    if logmdf:
        dn_dofe[dn_dofe == 0] = 1e-10
        weights = np.log10(dn_dofe)
    else:
        weights = dn_dofe
    axs[2].hist(ofe_bins[:-1], ofe_bins, weights=weights, color=color,
                orientation='horizontal', histtype=histtype, **style_kw)
    # ofe_bin_centers = get_bin_centers(ofe_bins)
    # axs[2].plot(weights, ofe_bin_centers, color=color, **style_kw)

    return fig, axs


def setup_axes(width=3.25, logmdf=False):
    """
    Create a figure with three axes: the main abundance track axis plus two
    side panels for [Fe/H] and [O/Fe] distribution functions.

    Parameters
    ----------
    width : float, optional
        Width of the figure in inches. The default is 3.25 in.
    logmdf : bool, optional
        If True, plot the marginal distributions on a log scale. The default
        is False.

    Returns
    -------
    fig : matplotlib.figure
    axs : list of matplotlib.axes.Axes
        Ordered [ax_main, ax_mdf, ax_odf]
    """
    fig = plt.figure(figsize=(width, width))
    gs = fig.add_gridspec(2, 2, width_ratios=(3, 1), height_ratios=(1, 3),
                          top=0.98, right=0.98, bottom=0.09, left=0.13,
                          wspace=0.05, hspace=0.05)
    # Start with the center panel for [Fe/H] vs [O/Fe]
    ax_main = fig.add_subplot(gs[1,0])
    ax_main.xaxis.set_major_locator(MultipleLocator(0.5))
    ax_main.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax_main.yaxis.set_major_locator(MultipleLocator(0.1))
    ax_main.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax_main.set_xlabel('[Fe/H]')
    ax_main.set_ylabel('[O/Fe]')
    ax_main.set_xlim((-2.5, 0.3))
    ax_main.set_ylim((-0.1, 0.54))
    # Add panel above for MDF in [Fe/H]
    ax_mdf = fig.add_subplot(gs[0,0], sharex=ax_main)
    ax_mdf.tick_params(axis='x', labelcolor='#ffffff00')
    ax_mdf.tick_params(axis='y', labelsize=7)
    if logmdf:
        ax_mdf.set_ylim((-4, 1.5))
        ax_mdf.yaxis.set_major_locator(MultipleLocator(2))
        ax_mdf.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax_mdf.set_ylabel(r'log(d$N$/d[Fe/H])', size=7)
    else:
        ax_mdf.set_ylabel(r'd$N$/d[Fe/H]', size=7)
        ax_mdf.yaxis.set_major_locator(MultipleLocator(5))
        ax_mdf.yaxis.set_minor_locator(MultipleLocator(1))
    # Add panel to the right for MDF in [O/Fe]
    ax_odf = fig.add_subplot(gs[1,1], sharey=ax_main)
    ax_odf.tick_params(axis='y', labelcolor='#ffffff00')
    ax_odf.tick_params(axis='x', labelsize=7)
    if logmdf:
        ax_odf.set_xlabel(r'log(d$N$/d[O/Fe])', size=7)
        ax_odf.set_xlim((-3.5, 1.5))
        ax_odf.xaxis.set_major_locator(MultipleLocator(2))
        ax_odf.xaxis.set_minor_locator(MultipleLocator(0.5))
    else:
        ax_odf.set_xlabel(r'd$N$/d[O/Fe]')
        ax_odf.xaxis.set_major_locator(MultipleLocator(5))
        ax_odf.xaxis.set_minor_locator(MultipleLocator(1))
    axs = [ax_main, ax_mdf, ax_odf]
    return fig, axs
