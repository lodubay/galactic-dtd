"""
This script plots abundance tracks for one-zone models with varying minimum
Type Ia delay times and star formation efficiency timescales.
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

# One-zone model settings
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    delay=0.04,
    tau_star=2.,
)

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'greggio05'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()


    simtime = np.arange(0, END_TIME + DT, DT)

    name = 'exponential30'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=3),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Exponential ($\tau=3$ Gyr)',
                      color=paultol.muted.colors[7],
                      style_kw={
                           'linestyle': '--',
                           'linewidth': 1,
                           'zorder': 1,
                      })

    name = 'exponential15'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=1.5),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Exponential ($\tau=1.5$ Gyr)',
                      color=paultol.muted.colors[6],
                      style_kw={
                           'linestyle': '--',
                           'linewidth': 1,
                           'zorder': 1,
                      })

    name = 'powerlaw'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.powerlaw(slope=-1.1),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Power-Law ($\alpha=-1.1$)',
                      color='k',
                      style_kw={
                           'linestyle': '--',
                           'linewidth': 1,
                           'zorder': 1,
                      })

    # Single-degenerate
    name = 'single'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.greggio05_single(),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label='Single degenerate',
                      color=paultol.muted.colors[0],
                       style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 9},
                      )

    # Double-degenerate WIDE model
    name = 'double_wide'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.greggio05_approximate.from_defaults('wide'),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label='Double degenerate WIDE',
                      color=paultol.muted.colors[1],
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 9},
                      )

    # Double-degenerate WIDE model
    name = 'double_close'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.greggio05_approximate.from_defaults('close'),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label='Double degenerate CLOSE',
                      color=paultol.muted.colors[3],
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 9},
                      )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.54))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_greggio05.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
