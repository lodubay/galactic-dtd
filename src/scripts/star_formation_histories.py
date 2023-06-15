"""
This script plots various quantities relating to star formation history as a
function of time for multiple VICE simulation outputs.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice
import paths
from utils import get_color_list, get_bin_centers
from _globals import END_TIME, ZONE_WIDTH, TWO_COLUMN_WIDTH, GALR_BINS
from colormaps import paultol
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.vibrant.colors)

MIGRATION = 'post-process'
EVOLUTION_LIST = ['insideout', 'lateburst', 'conroy22_JW20yields', 'twoinfall']
DTD = 'powerlaw_slope11'
CMAP_NAME = 'plasma_r'

def main():
    fig, axs = setup_axes()
    # Get color list
    cmap = plt.get_cmap(CMAP_NAME)
    colors = get_color_list(cmap, GALR_BINS)
    galr_centers = get_bin_centers(GALR_BINS)
    zones = [int(galr / ZONE_WIDTH) for galr in galr_centers]
    for i, evolution in enumerate(EVOLUTION_LIST):
        output = paths.data / 'migration' / MIGRATION / evolution / DTD
        for zone, color in zip(zones, colors):
            plot_history(axs[:,i], output, zone, color=color, 
                         label='%d kpc' % (zone * ZONE_WIDTH))
    leg = axs[0,0].legend(loc='upper center', frameon=False, ncols=3, 
                          handlelength=1, columnspacing=1, handletextpad=0.5)
    # Set legend text colors
    # for i in range(len(zones)):
    #     leg.get_texts()[i].set_color(colors[i])
    #     leg.legend_handles[i].set_visible(False)
    plt.savefig(paths.figures / 'star_formation_histories.pdf', dpi=300)
    plt.close()


def plot_history(axs, output, zone, color=None, label=None, linestyle='-',
                 zone_width=ZONE_WIDTH):
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
    radius = (zone + 0.5) * zone_width
    area = np.pi * ((radius + zone_width)**2 - radius**2)
    history = vice.history(str(Path('%s.vice' % output) / ('zone%i' % zone)))
    time = np.array(history['time'])
    sf_surface = np.array(history['sfr']) / area
    axs[0].plot(time, sf_surface, color=color, ls=linestyle, label=label)
    infall_surface = np.array(history['ifr']) / area
    axs[1].plot(time, infall_surface, color=color, label=label, ls=linestyle)
    gas_surface = np.array(history['mgas']) / area
    axs[2].plot(time, gas_surface, color=color, ls=linestyle, label=label)
    tau_star = [history['mgas'][i+1] / history['sfr'][i+1] * 1e-9 for i in range(
                len(history['time']) - 1)]
    axs[3].plot(history['time'][1:], tau_star, color=color, ls=linestyle,
                label=label)


def setup_axes(tmax=END_TIME, width=TWO_COLUMN_WIDTH):
    fig, axs = plt.subplots(4, 4, sharex=True, sharey='row', figsize=(width, width))
    fig.subplots_adjust(hspace=0., wspace=0., left=0.08, right=0.98, 
                        bottom=0.05, top=0.96)
    # Axis titles are model names
    axs[0,0].set_title('Inside-Out')
    axs[0,1].set_title('Late-Burst')
    axs[0,2].set_title('Early-Burst')
    axs[0,3].set_title('Two-Infall')
    # y-axis labels are diagnostics
    axs[0,0].set_ylabel(r'$\dot \Sigma_*$ [$M_{\odot}\,\rm{yr}^{-1}\,\rm{kpc}^{-2}$]')
    axs[1,0].set_ylabel(r'$\dot \Sigma_{\rm in}$ [$M_{\odot}\,\rm{yr}^{-1}\,\rm{kpc}^{-2}$]')
    axs[2,0].set_ylabel(r'$\Sigma_{\rm gas}$ [$M_{\odot}\,\rm{kpc}^{-2}$]')
    axs[3,0].set_ylabel(r'$\tau_*$ [Gyr]')
    # One x-axis label for all panels
    bigax = fig.add_subplot(111, frameon=False)
    bigax.tick_params(labelcolor='none', which='both', top=False, bottom=False, 
                      left=False, right=False)
    bigax.set_xlabel('Time [Gyr]')
    axs[0,0].set_xlim((-1, tmax+1))
    axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    # y-axis log scale
    axs[0,0].set_yscale('log')
    # axs[0,0].set_ylim((5e-5, 0.2))
    axs[0,0].set_ylim((5e-5, 0.5))
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylim((1e-3, 0.3))
    axs[2,0].set_yscale('log')
    axs[2,0].set_ylim((3e6, 3e8))
    axs[3,0].set_yscale('log')
    axs[3,0].set_ylim((0.5, 3e2))

    return fig, axs

if __name__ == '__main__':
    main()
