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
SLOPES = [-0.8, -1.1,-1.4]
TIMESCALES= [6, 3, 1.5]
DT = 0.01
DELAY = 0.04
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

# Plot settings
LINE_STYLE = ['-', '--', ':']
COLOR = ['k', 'k', 'k']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'slope_timescale'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=False)

    simtime = np.arange(0, END_TIME + DT, DT)

    for i, timescale in enumerate(TIMESCALES):
        dist = dtds.exponential(timescale=timescale, tmin=DELAY)
        # name = 'exponential{:02d}'.format(int(10*timescale))
        sz = vice.singlezone(name=str(output_dir / dist.name),
                             RIa=dist,
                             **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                          label=rf'Exponential ($\tau={timescale:.1f}$ Gyr)',
                          color=paultol.bright.colors[0],
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1},
                          marker_labels=(i==2),
                          logmdf=False
                          )

    for i, slope in enumerate(SLOPES):
        dist = dtds.powerlaw(slope=slope, tmin=DELAY)
        # name = 'powerlaw{:02d}'.format(int(-10*slope))
        sz = vice.singlezone(name=str(output_dir / dist.name),
                             RIa=dist,
                             **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                          label=rf'Power-Law ($\alpha={slope:.1f}$)',
                          color='k',
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1},
                          logmdf=False
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_slope_timescale.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
