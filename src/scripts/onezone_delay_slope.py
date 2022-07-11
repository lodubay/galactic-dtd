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
SLOPES = [-0.8, -1.1, -1.4]
TIMESCALES= [1.5, 3, 6]
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

# Plot settings
LINE_STYLE = [':', '--', '-']
COLOR = ['k', 'k', 'k']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'delay_slope'

    fig, axs = setup_axes()

    slope = -1.1
    delay = 0.15
    name = 'powerlaw{:02d}_delay{:03d}'.format(int(-10*slope), int(delay*1000))
    label = rf'$\alpha={slope:.2f}$, $t_D={int(delay*1000)}$ Myr'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.powerlaw(slope=slope, tmin=delay),
                         delay=delay,
                         **STANDARD_PARAMS)
    simtime = np.arange(0, END_TIME + DT, DT)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      plot_kw={'label': label},
                      style_kw={
                          'linestyle': '-',
                          'color': 'k',
                          'linewidth': 1,
                          'zorder': 10},
                      )

    timescale = 1.5
    delay = 0.04
    name = 'exponential{:02d}'.format(int(10*timescale))
    label = rf'Exponential ($\tau={timescale:.1f}$ Gyr)'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=timescale, tmin=delay),
                         delay=delay,
                         **STANDARD_PARAMS)
    simtime = np.arange(0, END_TIME + DT, DT)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      plot_kw={'label': label},
                      style_kw={
                          'linestyle': '-',
                          'color': paultol.bright.colors[5],
                          'linewidth': 1.5,
                          'zorder': 1},
                      )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.52))
    mdf_ylim = axs[1].get_ylim()
    axs[1].set_ylim((None, mdf_ylim[1]*2))
    odf_xlim = axs[2].get_xlim()
    axs[2].set_xlim((odf_xlim[0]*0.5, odf_xlim[1]*2))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_delay_slope.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
