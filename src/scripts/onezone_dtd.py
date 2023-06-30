"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
from multizone.src.simulations import models, dtds
from _globals import END_TIME
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
    eta=2.1,
    tau_star=2.,
    delay=DELAY,
)

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    distributions = [dtds.triple(tmin=DELAY),
                     dtds.exponential(timescale=3, tmin=DELAY),
                     dtds.plateau(width=0.3, slope=-1.1, tmin=DELAY),
                     dtds.powerlaw(slope=-1.1, tmin=DELAY),
                     dtds.prompt(peak=0.05, stdev=0.015, timescale=3, tmin=DELAY),
                     ]
    labels = [r'Triple-system evolution',
              r'Exponential ($\tau=3$ Gyr)',
              r'Plateau ($W=0.3$ Gyr)',
              r'Power law ($\alpha=-1.1$)',
              r'Prompt',]
    colors = [paultol.vibrant.colors[i] for i in [5, 1, 0, 4, 2]]
    line_styles = ['-', '-.', '--', ':', ':']

    for i, dist in enumerate(distributions):
        run_singlezone(str(output_dir / dist.name), simtime,
                       overwrite=overwrite,
                       RIa=dist, **STANDARD_PARAMS)
        plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                          label=labels[i], color=colors[i],
                          style_kw={'linestyle': '-',
                                    'linewidth': 1},
                          marker_labels=(i==0)
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.6))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd.png', dpi=300)
    fig.savefig(paths.figures / 'onezone_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
