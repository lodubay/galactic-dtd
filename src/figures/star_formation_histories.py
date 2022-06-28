"""
This script plots various quantities relating to star formation history as a
function of time for multiple VICE simulation outputs.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice

def main():
    fig, axs = setup_axes()
    plot_history('../data/migration_outputs/post-process/insideout/powerlaw', 80, axs)
    plot_history('../data/migration_outputs/post-process/conroy22_insideout/powerlaw', 80, axs)
    plt.show()

def plot_history(output, zone, axs):
    r"""
    Plot IFR, SFR, Mgas, and SFE timescale for the given VICE multioutput and
    zone.

    Parameters
    ----------
    output : str
        Path to multioutput directory.
    zone : int
        Index of zone to plot.
    axs : list of Axes

    """
    history = vice.history(str(Path(output + '.vice') / ('zone%i' % zone)))
    axs[0,0].plot(history['time'], history['ifr'])
    axs[0,1].plot(history['time'], history['sfr'])
    axs[1,0].plot(history['time'], history['mgas'])
    tau_star = [history['mgas'][i+1] / history['sfr'][i+1] * 1e-9 for i in range(
                len(history['time']) - 1)]
    axs[1,1].plot(history['time'][1:], tau_star)


def setup_axes():
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(8, 8))

    axs[0,0].set_title(r'Infall Rate [$M_{\odot}\,\rm{yr}^{-1}$]')
    axs[0,0].set_xlim((0, 13.2))
    axs[0,1].set_title(r'Star Formation Rate [$M_{\odot}\,\rm{yr}^{-1}$]')
    axs[1,0].set_title(r'Gas Mass [$M_{\odot}$]')
    axs[1,1].set_title(r'SFE Timescale [Gyr]')

    return fig, axs

if __name__ == '__main__':
    main()