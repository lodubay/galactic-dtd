import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def main():
    import vice
    import paths
    import sys
    sys.path.append(str(paths.root))
    from migration.src.simulations import models, dtds
    from vice.yields.presets import JW20

    fname = str(paths.data / 'onezone' / 'test01')
    sz = vice.singlezone(name=fname, func=models.insideout(8), mode='sfr',
                         RIa=dtds.powerlaw(), delay=0.04)
    simtime = np.arange(0, 13.21, 0.01)
    sz.run(simtime, overwrite=True)
    hist = vice.history(fname)
    mdf = vice.mdf(fname)
    mdf_bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
    fig, axs = plot_track_and_mdf(hist['[fe/h]'], hist['[o/fe]'],
                                  dn_dfeh=mdf['dn/d[fe/h]'], feh_bins=mdf_bins,
                                  dn_dofe=mdf['dn/d[o/fe]'], ofe_bins=mdf_bins,
                                  feh_lim=(-2.8, 0.), ofe_lim=(-0.1, 0.5),
                                  filled=True)

    fname = str(paths.data / 'onezone' / 'test02')
    sz = vice.singlezone(name=fname, func=models.insideout(8), mode='sfr',
                         RIa=dtds.bimodal(), delay=0.04)
    simtime = np.arange(0, 13.21, 0.01)
    sz.run(simtime, overwrite=True)
    hist = vice.history(fname)
    mdf = vice.mdf(fname)
    mdf_bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
    plot_track_and_mdf(hist['[fe/h]'], hist['[o/fe]'],
                       dn_dfeh=mdf['dn/d[fe/h]'], feh_bins=mdf_bins,
                       dn_dofe=mdf['dn/d[o/fe]'], ofe_bins=mdf_bins,
                       fig=fig, axs=axs,
                       feh_lim=(-2.8, 0.), ofe_lim=(-0.1, 0.5),
                       filled=False, color='r')

    plt.savefig(paths.figures / 'track_dist_test.png', dpi=300)
    plt.close()


def plot_vice_onezone(output, fig=None, axs=[], **kwargs):
    import vice
    hist = vice.history(output)
    mdf = vice.mdf(output)
    mdf_bins = mdf['bin_edge_left'] + mdf['bin_edge_right'][-1:]
    fig, axs = plot_track_and_mdf(hist['[fe/h]'], hist['[o/fe]'],
                       dn_dfeh=mdf['dn/d[fe/h]'], feh_bins=mdf_bins,
                       dn_dofe=mdf['dn/d[o/fe]'], ofe_bins=mdf_bins,
                       fig=fig, axs=axs, **kwargs)
    return fig, axs



def plot_track_and_mdf(feh, ofe, dn_dfeh=[], feh_bins=[], dn_dofe=[], ofe_bins=[],
                       fig=None, axs=[], color='k', linestyle='-',
                       linewidth=1, zorder=10, filled=False, label=None,
                       **kwargs):
    if fig == None:
        fig, axs = setup_axes(**kwargs)
    if filled:
        histtype = 'stepfilled'
    else:
        histtype = 'step'
    axs[0].plot(feh, ofe, color=color, linestyle=linestyle, linewidth=linewidth,
                zorder=zorder, label=label)
    axs[1].hist(feh_bins[:-1], feh_bins, weights=dn_dfeh,
                histtype=histtype, color=color, linestyle=linestyle,
                linewidth=linewidth, zorder=zorder)
    axs[2].hist(ofe_bins[:-1], ofe_bins, weights=dn_dofe,
                orientation='horizontal',
                histtype=histtype, color=color, linestyle=linestyle,
                linewidth=linewidth, zorder=zorder)
    return fig, axs


def setup_axes(feh_lim=(-3, 1), ofe_lim=(-0.2, 0.6), width=3.25, **kwargs):
    """
    Create a figure with three axes: the main abundance track axis plus two
    side panels for [Fe/H] and [O/Fe] distribution functions.

    Parameters
    ----------
    feh_lim : tuple
        Bounds of the x- or [Fe/H]-axis. The default is (-3, 1).
    ofe_lim : tuple
        Bounds of the y- or [O/Fe]-axis. The default is (-0.2, 0.6).
    width : float
        Width of the figure in inches. The default is 3.25 in.
    Other keyword arguments are passed to plt.subplots

    Returns
    -------
    fig : matplotlib.figure
    axs : list of matplotlib.axes.Axes
    """
    # Start with the center panel
    fig, ax_main = plt.subplots(figsize=(width, width), **kwargs)
    # fig.subplots_adjust(left=0.3, bottom=0.25, right=0.95, top=0.95)
    fig.subplots_adjust(left=0.15, bottom=0.12, right=0.9, top=0.95)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.tick_params(top=False, right=False, which='both')
    ax_main.tick_params(direction='out', which='both')
    ax_main.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax_main.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax_main.set_xlim(feh_lim)
    ax_main.set_ylim(ofe_lim)
    ax_main.set_xlabel('[Fe/H]')
    ax_main.set_ylabel('[O/Fe]')
    # Add panel below for MDF
    # ax_mdf = plt.axes([0.3, 0.1, 0.65, 0.12], sharex=ax_main)
    ax_mdf = plt.axes([0.15, 0.12, 0.75, 0.1], sharex=ax_main)
    ax_mdf.set_yscale('log')
    # ax_mdf.spines['top'].set_visible(False)
    # ax_mdf.spines['bottom'].set_visible(False)
    # ax_mdf.spines['left'].set_visible(False)
    # ax_mdf.spines['right'].set_visible(False)
    # ax_mdf.tick_params(top=False, left=False, bottom=False, right=True,
    #                    which='both', direction='out')
    # ax_mdf.minorticks_off()
    # ax_mdf.tick_params(axis='x', labelcolor='#ffffff00')
    # ax_mdf.tick_params(axis='y', labelsize=6)
    # ax_mdf.yaxis.set_label_position('right')
    # ax_mdf.yaxis.tick_right()
    # ax_mdf.set_ylim((0.0001, 10))
    ax_mdf.axis('off')
    # Add panel to the left for MDF in [O/Fe]
    # ax_odf = plt.axes([0.15, 0.25, 0.12, 0.7], sharey=ax_main)
    ax_odf = plt.axes([0.15, 0.12, 0.1, 0.83], sharey=ax_main)
    ax_odf.set_xscale('log')
    ax_odf.axis('off')
    axs = [ax_main, ax_mdf, ax_odf]
    return fig, axs


if __name__ == '__main__':
    main()

