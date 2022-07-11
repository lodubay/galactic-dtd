"""
This file contains functions for plotting abundance tracks in [Fe/H] vs [O/Fe]
alongside their corresponding metallicity distribution functions (MDFs).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, LogLocator

def main():
    """
    Proof-of-concept plot
    """
    import vice
    import paths
    import sys
    sys.path.append(str(paths.root))
    from migration.src.simulations import models, dtds
    from vice.yields.presets import JW20

    fig, axs = setup_axes()

    fname = str(paths.data / 'onezone' / 'test01')
    sz = vice.singlezone(name=fname, func=models.insideout(8), mode='sfr',
                         RIa=dtds.powerlaw(), delay=0.04)
    simtime = np.arange(0, 13.21, 0.01)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(fname, fig=fig, axs=axs, style_kw={'color': 'k'},
                      hist_kw={'histtype': 'stepfilled'})

    fname = str(paths.data / 'onezone' / 'test02')
    sz = vice.singlezone(name=fname, func=models.insideout(8), mode='sfr',
                         RIa=dtds.exponential(), delay=0.04)
    simtime = np.arange(0, 13.21, 0.01)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(fname, fig=fig, axs=axs, hist_kw={'histtype': 'step'},
                      style_kw={'color': 'r'})

    axs[0].set_xlim((-2.6, 0.1))
    axs[0].set_ylim((-0.02, 0.52))

    plt.savefig(paths.figures / 'track_dist_test.png', dpi=300)
    plt.close()


def plot_vice_onezone(output, **kwargs):
    """
    Wrapper for plot_track_and_mdf given a VICE onezone output.

    Parameters
    ----------
    output : str
        Path to VICE output, not necessarily including the '.vice' extension
    **kwargs passed to plot_track_and_mdf

    Returns
    -------
    fig, ax
    """
    import vice
    hist = vice.history(output)
    mdf = vice.mdf(output)
    mdf_bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
    fig, axs = plot_track_and_mdf(hist['[fe/h]'], hist['[o/fe]'],
                                  dn_dfeh=mdf['dn/d[fe/h]'], feh_bins=mdf_bins,
                                  dn_dofe=mdf['dn/d[o/fe]'], ofe_bins=mdf_bins,
                                  **kwargs)
    return fig, axs


def plot_track_and_mdf(feh, ofe, dn_dfeh=[], feh_bins=None, dn_dofe=[],
                       ofe_bins=None, fig=None, axs=[], style_kw={},
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
    feh_bins : array-like, optional
        Bins in [Fe/H] for plotting the MDF in [Fe/H]
    dn_dofe : array-like, optional
        Distribution of O abundances dN/d[O/Fe]. If empty, the ofe parameter
        is passed to matplotlib.pyplot.hist directly
    ofe_bins : array-like, optional
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
        Dict of style-related keyword arguments to pass to matplotlib.pyplot.plot
    hist_kw : dict, optional
        Dict of style-related keyword arguments to pass to matplotlib.pyplot.plot

    Returns
    -------
    fig, axs
    """
    if fig == None or len(axs) != 3:
        fig, axs = setup_axes()

    axs[0].plot(feh, ofe, **plot_kw, **style_kw)
    if len(dn_dfeh) > 0:
        axs[1].hist(feh_bins[:-1], feh_bins, weights=dn_dfeh,
                    **hist_kw, **style_kw)
    else:
        axs[1].hist(feh, bins=feh_bins, density=True, **hist_kw, **style_kw)
    if len(dn_dofe) > 0:
        axs[2].hist(ofe_bins[:-1], ofe_bins, weights=dn_dofe,
                    orientation='horizontal', **hist_kw, **style_kw)
    else:
        axs[2].hist(ofe, bins=ofe_bins, density=True, orientation='horizontal',
                    **hist_kw, **style_kw)
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
    ax_mdf.set_yscale('log')
    ax_mdf.set_ylabel('dN/d[Fe/H]')
    ax_mdf.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
    # Add panel to the right for MDF in [O/Fe]
    ax_odf = fig.add_subplot(gs[1,1], sharey=ax_main)
    ax_odf.tick_params(axis='y', labelcolor='#ffffff00')
    ax_odf.set_xscale('log')
    ax_odf.set_xlabel('dN/d[O/Fe]')
    ax_odf.xaxis.set_major_locator(LogLocator(base=10, numticks=5))
    axs = [ax_main, ax_mdf, ax_odf]
    return fig, axs


if __name__ == '__main__':
    main()

