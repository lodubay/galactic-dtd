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

    # Double-degenerate WIDE model
    dtd = dtds.greggio05_approximate.from_defaults('wide')
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Double Degenerate WIDE',
                      color=paultol.muted.colors[1],
                      marker_labels=True,
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 10},
                      )
    
    # Double-degenerate CLOSE model
    dtd = dtds.greggio05_approximate.from_defaults('close')
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Double Degenerate CLOSE',
                      color=paultol.muted.colors[3],
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 9},
                      )

    # Single-degenerate
    dtd = dtds.greggio05_single()
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Single Degenerate',
                      color=paultol.muted.colors[0],
                       style_kw={
                           'linestyle': '-',
                           'linewidth': 1,
                           'zorder': 8},
                      )

    dtd = dtds.exponential(timescale=3)
    sz = vice.singlezone(name=str(output_dir / dtd.name), 
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Exponential ($\tau=3$ Gyr)',
                      color=paultol.muted.colors[7],
                      style_kw={
                            'linestyle': '--',
                            'linewidth': 1,
                            'zorder': 1,
                      })

    dtd = dtds.exponential(timescale=1.5)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Exponential ($\tau=1.5$ Gyr)',
                      color=paultol.muted.colors[6],
                      style_kw={
                           'linestyle': '--',
                           'linewidth': 1,
                           'zorder': 1,
                      })

    dtd = dtds.powerlaw(slope=-1.1)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Power-Law ($\alpha=-1.1$)',
                      color='k',
                      style_kw={
                           'linestyle': '--',
                           'linewidth': 1,
                           'zorder': 1,
                      })

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.54))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_greggio05.pdf', dpi=300)
    plt.close()    


if __name__ == '__main__':
    main()
