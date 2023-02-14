"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import models, dtds
from migration.src._globals import END_TIME
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from utils import run_singlezone

# VICE one-zone model settings
DELAYS = [0.04, 0.15]
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    tau_star=2.,
)
LINE_STYLES = ['-', '--']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'dtd_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    labels = [r'Power-Law ($\alpha=-1.1$)',
              r'Power-Law with 200 Myr plateau',
              r'Exponential ($\tau=3$ Gyr)',]
    colors = [paultol.vibrant.colors[i] for i in [4, 0, 1, 2]]

    for delay, ls in zip(DELAYS, LINE_STYLES):
        distributions = [dtds.powerlaw(slope=-1.1, tmin=delay),
                         dtds.plateau(width=0.3, slope=-1.1, tmin=delay),
                         dtds.exponential(timescale=3, tmin=delay)]
        for i, dist in enumerate(distributions):
            if delay == DELAYS[0]:
                label = labels[i]
            else:
                label = None
            name = dist.name + f'_delay{int(1000*delay)}'
            run_singlezone(str(output_dir / name), simtime,
                           overwrite=overwrite,
                           RIa=dist, delay=delay, **STANDARD_PARAMS)
            plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                              label=label, color=colors[i],
                              style_kw={'linestyle': ls,
                                        'linewidth': 1},
                              marker_labels=(i==2 and delay==DELAYS[1])
                              )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.52))

    # Legend
    handles, labels = axs[0].get_legend_handles_labels()
    handles = [Line2D([], [], color='k', ls=ls, lw=1) for ls in LINE_STYLES] +\
        handles
    labels = [f'{int(1000*delay)} Myr minimum delay' for delay in DELAYS] + \
        labels
    axs[0].legend(handles, labels, frameon=False, loc='lower left',
                  handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd_delay.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
