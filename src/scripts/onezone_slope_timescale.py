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

# VICE one-zone model settings
SLOPES = [-1.4, -1.1, -0.8]
TIMESCALES= [6, 3, 1.5]
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
LINE_STYLE = [':', '--', '-']
COLOR = ['k', 'k', 'k']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'slope_timescale'

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    for i, timescale in enumerate(TIMESCALES):
        name = 'exponential{:02d}'.format(int(10*timescale))
        label = rf'Exponential ($\tau={timescale:.1f}$ Gyr)'
        sz = vice.singlezone(name=str(output_dir / name),
                             RIa=dtds.exponential(timescale=timescale),
                             **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                          plot_kw={'label': label},
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'color': paultol.bright.colors[5],
                              'linewidth': 1},
                          )

    for i, slope in enumerate(SLOPES):
        name = 'powerlaw{:02d}'.format(int(-10*slope))
        label = rf'Power-Law ($\alpha={slope:.1f}$)'
        sz = vice.singlezone(name=str(output_dir / name),
                             RIa=dtds.powerlaw(slope=slope),
                             **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                          plot_kw={'label': label},
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'color': 'k',
                              'linewidth': 1},
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_slope_timescale.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
