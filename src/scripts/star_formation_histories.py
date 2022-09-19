"""
This script plots various quantities relating to star formation history as a
function of time for multiple VICE simulation outputs.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from utils import discrete_colormap, get_color_list, setup_discrete_colorbar, \
    get_bin_centers
from _globals import END_TIME, GALR_BINS, ZONE_WIDTH

def main(evolution, RIa='powerlaw_slope11', cmap_name='plasma_r'):
    output = paths.data / 'migration' / 'post-process' / evolution / RIa
    fname = 'sfh_%s.png' % evolution
    plot_by_zone(output, fname=fname, zone_width=ZONE_WIDTH, cmap_name=cmap_name)
    # fig, axs = setup_axes()
    # plot_history(axs, paths.data / 'migration/diffusion/twoinfall/powerlaw_slope11_delay040', 80, color='b')
    # plot_history(axs, paths.data / 'migration/diffusion/insideout_conroy22/powerlaw_slope11_delay040', 80, color='orange')
    # plt.savefig(paths.figures / 'star_formation_histories.png', dpi=300)
    

def plot_by_zone(output, galr_bins=GALR_BINS, zone_width=ZONE_WIDTH, 
                 cmap_name='plasma_r', fname='star_formation_histories.png'):
    # Determine zones to plot
    galr_centers = get_bin_centers(galr_bins)
    zones = [int(galr / zone_width) for galr in galr_centers]
    # Set up plot
    fig, axs = setup_axes()
    cmap, norm = discrete_colormap(cmap_name, galr_bins)
    colors = get_color_list(cmap, galr_bins)
    # Plot SFHs for each zone
    for zone, color in zip(zones, colors):
        plot_history(axs, output, zone, color=color, 
                     label='%s kpc' % int(zone * zone_width))
    axs[0,0].legend(frameon=False, fontsize=8, title='Galactic radius', 
                    title_fontsize=8)
    plt.savefig(paths.figures / fname, dpi=300)
    plt.close()


def plot_history(axs, output, zone, color=None, label=None):
    r"""
    Plot IFR, SFR, Mgas, and SFE timescale for the given VICE multioutput and
    zone.

    Parameters
    ----------
    axs : list of Axes
    output : str
        Path to multioutput directory.
    zone : int
        Index of zone to plot.
    color : str or None, optional
        Plot color. The default is None.
    """
    history = vice.history(str(Path('%s.vice' % output) / ('zone%i' % zone)))
    axs[0,0].plot(history['time'], history['ifr'],
                  color=color, label=label)
    axs[0,1].plot(history['time'], history['sfr'], color=color)
    axs[1,0].plot(history['time'], history['mgas'], color=color)
    tau_star = [history['mgas'][i+1] / history['sfr'][i+1] * 1e-9 for i in range(
                len(history['time']) - 1)]
    axs[1,1].plot(history['time'][1:], tau_star, color=color)


def setup_axes(tmax=END_TIME):
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(7, 7))

    axs[0,0].set_title(r'Infall Rate [$M_{\odot}\,\rm{yr}^{-1}$]')
    axs[0,0].set_xlim((-1, tmax+1))
    axs[0,1].set_title(r'Star Formation Rate [$M_{\odot}\,\rm{yr}^{-1}$]')
    axs[1,0].set_title(r'Gas Mass [$M_{\odot}$]')
    axs[1,0].set_xlabel('Time [Gyr]')
    axs[1,1].set_title(r'Star Formation Efficiency Timescale [Gyr]')
    axs[1,1].set_yscale('log')
    axs[1,1].set_xlabel('Time [Gyr]')

    return fig, axs

if __name__ == '__main__':
    evolution = sys.argv[1]
    main(evolution)
