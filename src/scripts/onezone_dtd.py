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
from migration.src.simulations import models
from migration.src._globals import END_TIME
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

# VICE one-zone model settings
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    tau_star=2.,
    delay=0.04,
)

# Plot settings
DTDS = [styles.bimodal, styles.plaw_steep, styles.plaw, styles.plaw_broken,
        styles.exp, styles.exp_long]

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'dtd'

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    for dtd in DTDS:
        sz = vice.singlezone(name=str(output_dir / dtd['name']),
                             RIa=dtd['func'], **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd['name']), fig=fig, axs=axs,
                          plot_kw={'label': dtd['label']},
                          style_kw={
                              'linestyle': dtd['line'],
                              'color': dtd['color'],
                              'linewidth': 1},
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
