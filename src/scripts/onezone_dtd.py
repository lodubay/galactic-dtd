"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
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
DELAY = 0.04
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    tau_star=2.,
    delay=DELAY,
)

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    distributions = [dtds.powerlaw(slope=-1.1, tmin=DELAY),
                     dtds.plateau(width=0.2, slope=-1.1, tmin=DELAY),
                     dtds.exponential(timescale=3, tmin=DELAY),
                     dtds.prompt(center=0.05, timescale=3, tmin=DELAY)]
    labels = [r'Power-Law ($\alpha=-1.1$)',
              r'Power-Law with 200 Myr plateau',
              r'Exponential ($\tau=3$ Gyr)',
              r'Exponential with prompt component',]
    colors = [paultol.vibrant.colors[i] for i in [4, 0, 1, 2]]
    line_styles = ['-', '-.', '--', ':']

    for i, dist in enumerate(distributions):
        run_singlezone(str(output_dir / dist.name), simtime,
                       overwrite=overwrite,
                       RIa=dist, **STANDARD_PARAMS)
        plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                          label=labels[i], color=colors[i],
                          style_kw={'linestyle': line_styles[i],
                                    'linewidth': 1},
                          marker_labels=(i==2)
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
