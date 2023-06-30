"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
LOG_MDF = False
SLOPE = -1.1

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'powerlaw_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=LOG_MDF)

    simtime = np.arange(0, END_TIME + DT, DT)

    delay = 0.04
    dist = dtds.powerlaw(slope=SLOPE, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Power-Law ($t_D=40$ Myr)', 
                      color='k',
                      style_kw={'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 10},
                      logmdf=LOG_MDF
                      )
    
    delay = 0.15
    dist = dtds.powerlaw(slope=SLOPE, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Power-Law ($t_D=150$ Myr)', 
                      color='k',
                      style_kw={'linestyle': '--',
                                'linewidth': 1,
                                'zorder': 10},
                      marker_labels=True,
                      logmdf=LOG_MDF
                      )
    
    delay = 0.04
    dist = dtds.plateau(width=0.35, slope=SLOPE, tmin=delay)
    name = dist.name + f'_delay{int(1000*delay)}'
    run_singlezone(str(output_dir / name), simtime,
                   overwrite=overwrite,
                   RIa=dist, delay=delay, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                      label=r'Plateau ($W=350$ Myr, $t_D=40$ Myr)', 
                      color='b',
                      style_kw={'linestyle': '-',
                                'linewidth': 1.5,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )
    
    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_powerlaw_delay.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
