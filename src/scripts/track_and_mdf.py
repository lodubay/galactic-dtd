"""
This file contains functions for plotting abundance tracks in [Fe/H] vs [O/Fe]
alongside their corresponding metallicity distribution functions (MDFs).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def plot_vice_onezone(output, fig=None, axs=[], style_kw={}, plot_kw={},
                      hist_kw={'histtype': 'step'}):
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
    style_kw : dict, optional
        Dict of style-related keyword arguments to pass to both
        matplotlib.pyplot.plot and matplotlib.pyplot.hist
    plot_kw : dict, optional
        Dict of style-related keyword arguments to pass to
        matplotlib.pyplot.plot
    hist_kw : dict, optional
        Dict of style-related keyword arguments to pass to
        matplotlib.pyplot.plot

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
                                  fig=fig, axs=axs, style_kw=style_kw,
                                  plot_kw=plot_kw, hist_kw=hist_kw)
    return fig, axs


def plot_track_and_mdf(feh, ofe, dn_dfeh=[], feh_bins=10, dn_dofe=[],
                       ofe_bins=10, fig=None, axs=[], style_kw={},
                       plot_kw={}, hist_kw={'histtype': 'step'}):
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
    style_kw : dict, optional
        Dict of style-related keyword arguments to pass to both
        matplotlib.pyplot.plot and matplotlib.pyplot.hist
    plot_kw : dict, optional
        Dict of style-related keyword arguments to pass to
        matplotlib.pyplot.plot
    hist_kw : dict, optional
        Dict of style-related keyword arguments to pass to
        matplotlib.pyplot.plot

    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : list of matplotlib.axes.Axes
    """
    if fig == None or len(axs) != 3:
        fig, axs = setup_axes()

    # Plot abundance tracks on main panel
    axs[0].plot(feh, ofe, **plot_kw, **style_kw)

    # Plot distribution of [Fe/H] on top panel
    if len(dn_dfeh) == 0:
        dn_dfeh, feh_bins = np.histogram(feh, bins=feh_bins)
    # mask zeros before taking log
    dn_dfeh = np.array(dn_dfeh)
    dn_dfeh[dn_dfeh == 0] = 1e-10
    log_dn_dfeh = np.log10(dn_dfeh)
    axs[1].hist(feh_bins[:-1], feh_bins, weights=log_dn_dfeh,
                **hist_kw, **style_kw)

    # Plot distribution of [O/Fe] on right side panel
    if len(dn_dofe) == 0:
        dn_dofe, ofe_bins = np.histogram(ofe, bins=ofe_bins)
    dn_dofe = np.array(dn_dofe)
    dn_dofe[dn_dofe == 0] = 1e-10
    log_dn_dofe = np.log10(dn_dofe)
    axs[2].hist(ofe_bins[:-1], ofe_bins, weights=log_dn_dofe,
                orientation='horizontal', **hist_kw, **style_kw)

    return fig, axs


def setup_axes(width=3.25):
    """
    Create a figure with three axes: the main abundance track axis plus two
    side panels for [Fe/H] and [O/Fe] distribution functions.

    Parameters
    ----------
    width : float
        Width of the figure in inches. The default is 3.25 in.

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
    ax_main.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax_main.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax_main.set_xlabel('[Fe/H]')
    ax_main.set_ylabel('[O/Fe]')
    # Add panel above for MDF in [Fe/H]
    ax_mdf = fig.add_subplot(gs[0,0], sharex=ax_main)
    ax_mdf.tick_params(axis='x', labelcolor='#ffffff00')
    ax_mdf.set_ylim((-4, 1.5))
    ax_mdf.yaxis.set_major_locator(MultipleLocator(2))
    ax_mdf.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_mdf.set_ylabel('log(dN/d[Fe/H])', size=7)
    ax_mdf.tick_params(axis='y', labelsize=7)
    # Add panel to the right for MDF in [O/Fe]
    ax_odf = fig.add_subplot(gs[1,1], sharey=ax_main)
    ax_odf.tick_params(axis='y', labelcolor='#ffffff00')
    # ax_odf.set_xscale('log')
    ax_odf.set_xlabel('log(dN/d[O/Fe])', size=7)
    ax_odf.tick_params(axis='x', labelsize=7)
    ax_odf.set_xlim((-3.5, 1.5))
    ax_odf.xaxis.set_major_locator(MultipleLocator(2))
    ax_odf.xaxis.set_minor_locator(MultipleLocator(0.5))
    axs = [ax_main, ax_mdf, ax_odf]
    return fig, axs
