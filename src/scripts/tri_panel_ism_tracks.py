"""
Plot tracks in [Fe/H] vs age, [O/Fe] vs age, and [O/Fe] vs [Fe/H] for the ISM
in the 8 kpc annulus from VICE multizone simulations for different DTDs
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice

standard = {
    'name': 'powerlaw',
    'label': r'$\alpha=-1.1$, $\tau_{\rm{min}}=40$ Myr',
    'color': 'k',
    'line': '-',
}
steep = {
    'name': 'powerlaw_steep',
    'label': r'$\alpha=-1.4$, $\tau_{\rm{min}}=40$ Myr',
    'color': '#4477aa',
    'line': '--',
}
standard_delayed = {
    'name': 'long_delay',
    'label': r'$\alpha=-1.1$, $\tau_{\rm{min}}=150$ Myr',
    'color': '#ee6677',
    'line': '-.',
}
steep_delayed = {
    'name': 'powerlaw_steep_delayed',
    'label': r'$\alpha=-1.4$, $\tau_{\rm{min}}=150$ Myr',
    'color': '#ccbb44',
    'line': ':',
}
exp = {
    'name': 'exponential',
    'label': r'Exponential ($\tau=1.5$ Gyr)',
    'color': '#228833',
    'line': '-',
}

def main():
    fig, axs = setup_axes()
    dtds = [standard, steep, standard_delayed, steep_delayed, exp]
    for dtd in dtds:
        plot_tracks(axs, dtd)
    axs[0].legend(fontsize=10, frameon=False)
    fig.suptitle('Power-Law DTD Variants')
    plt.savefig('tri_panel_ism_tracks.pdf', dpi=300)
    plt.close()


def plot_tracks(axs, dtd, zone=80,
                parent_dir='../data/migration_outputs/post-process/insideout'):
    """
    Plot abundance tracks from a particular zone (default: solar) in a VICE
    multizone simulation.

    Parameters
    ----------
    axs : list of axes
    dtd : dict
        One of the above delay-time distribution dictionaries
    zone : int, optional [default: 80]
        Zone number to plot
    parent_dir : str, optional
        Parent directory of VICE migration outputs
    """
    zone_path = Path(parent_dir) / Path('%s.vice/zone%s' % (dtd['name'], zone))
    hist = vice.history(str(zone_path))
    axs[0].plot(hist['time'], hist['[fe/h]'], color=dtd['color'],
                linestyle=dtd['line'], label=dtd['label'])
    axs[1].plot(hist['time'], hist['[o/fe]'], color=dtd['color'],
                linestyle=dtd['line'])
    axs[2].plot(hist['[fe/h]'], hist['[o/fe]'], color=dtd['color'],
                linestyle=dtd['line'])


def setup_axes(tlim=(-1, 14), felim=(-3, 0.25), olim=(-0.1, 0.5)):
    """
    Set up three-panel plot.

    Parameters
    ----------
    tlim : tuple, optional [default: (-1, 14)]
        Bounds of time axes
    felim : tuple, optional [default: (-2, 0.5)]
        Bounds of [Fe/H] axes
    olim : tuple, optional [default: (-0.1, 0.45)]
        Bounds of [O/Fe] axes

    Returns
    -------
    fig : figure
    axs : list of axes
    """
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    # First panel: metallicity vs time
    ax = axs[0]
    ax.set_xlim(tlim)
    ax.set_ylim(felim)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel('[Fe/H]')

    # Second panel: alpha vs time
    ax = axs[1]
    ax.set_xlim(tlim)
    ax.set_ylim(olim)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel('[O/Fe]')

    # Third panel: alpha vs fe
    ax = axs[2]
    ax.set_xlim(felim)
    ax.set_ylim(olim)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('[O/Fe]')

    return fig, axs

if __name__ == '__main__':
    main()
