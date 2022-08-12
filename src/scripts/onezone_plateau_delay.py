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
DELAYS = [0.15, 0.04]
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
LINE_STYLES = ['--', '-']
COLORS = paultol.highcontrast.colors
LOG_MDF = False

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'plateau_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=LOG_MDF)

    simtime = np.arange(0, END_TIME + DT, DT)

    for delay, ls in zip(DELAYS, LINE_STYLES):
        for width, color in zip()
    delay = 0.15
    dist = dtds.plateau(width=1, slope=-1.1, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=None, 
                      color='b',
                      style_kw={'linestyle': '--',
                                'linewidth': 1,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )
    
    delay = 0.04
    dist = dtds.plateau(width=1, slope=-1.1, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'1 Gyr plateau', 
                      color='b',
                      style_kw={'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )
                      
    delay = 0.04
    dist = dtds.plateau(width=0.35, slope=-1.1, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'350 Myr plateau', 
                      color='r',
                      style_kw={'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )
    
    delay = 0.15
    dist = dtds.plateau(width=0.1, slope=-1.1, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=None, 
                      color='orange',
                      style_kw={'linestyle': '--',
                                'linewidth': 1,
                                'zorder': 2},
                      logmdf=LOG_MDF
                      )
                  
    delay = 0.04
    dist = dtds.plateau(width=0.1, slope=-1.1, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'100 Myr plateau', 
                      color='orange',
                      style_kw={'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 2},
                      logmdf=LOG_MDF
                      )

    # Legend
    handles, labels = axs[0].get_legend_handles_labels()
    handles += [Line2D([], [], color='k', ls=ls, lw=1) for ls in LINE_STYLES]
    labels += [f'{int(1000*delay)} Myr minimum delay' for delay in DELAYS]
    axs[0].legend(handles, labels, frameon=False, loc='lower left',
                  handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_plateau_delay.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
