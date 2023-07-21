"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS, DT
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from utils import run_singlezone

# Set color cycle from Paul Tol colorscheme
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)

# DTD parameters
PEAKS = [0.5, 0.2, 0.1, 0.05]
# PEAKS = [0.5, 0.5, 0.5, 0.5]
# WIDTHS = [0.1, 0.04, 0.02, 0.01]
WIDTHS = [0.15, 0.06, 0.03, 0.015]
TIMESCALE = 3 # Exponential timescale in Gyr

# Plot settings
LINE_STYLE = [':', '-.', '--', '-']

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')

    output_dir = paths.data / 'onezone' / 'prompt'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=False)

    simtime = np.arange(0, END_TIME + DT, DT)

    # Exponential DTD for comparison
    dtd = dtds.exponential(timescale=TIMESCALE)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, 
                         func=models.insideout(8, dt=DT), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=rf'Exponential ($\tau={TIMESCALE:.01f}$ Gyr)',
                      style_kw={'linewidth': 1, 'zorder': 1},
                      marker_labels=True,
                      logmdf=False
    )

    for i in range(len(PEAKS)):
        dtd = dtds.prompt(peak=PEAKS[i], stdev=WIDTHS[i],
                          timescale=TIMESCALE)
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd, 
                             func=models.insideout(8, dt=DT), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                          label=r'$\bar t={}$ Myr, $\sigma_t={}$ Myr'.format(
                              int(PEAKS[i]*1000), int(WIDTHS[i]*1000)),
                          style_kw={'linestyle': LINE_STYLE[i], 
                                    'linewidth': 1, 
                                    'zorder': 10-i},
                          logmdf=False
        )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.4))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_twopopulation.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
