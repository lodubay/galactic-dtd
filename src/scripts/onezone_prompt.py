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
from utils import run_singlezone

# Set color cycle from Paul Tol colorscheme
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)

# VICE one-zone model settings
PEAKS = [0.5, 0.2, 0.1, 0.05]
WIDTHS = [0.1, 0.04, 0.02, 0.01]
TIMESCALE = 3 # Exponential timescale in Gyr
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
LINE_STYLE = [':', '-.', '--', '-']

def main(overwrite=False):

    output_dir = paths.data / 'onezone' / 'prompt'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    # Exponential DTD for comparison
    dist = dtds.exponential(timescale=TIMESCALE)
    run_singlezone(str(output_dir / dist.name), simtime, overwrite=overwrite,
                   RIa=dist, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                      label=rf'Exponential ($\tau={TIMESCALE:.01f}$ Gyr)',
    )

    for i in range(len(PEAKS)):
        dist = dtds.prompt(center=PEAKS[i], stdev=WIDTHS[i],
                           timescale=TIMESCALE, tmin=DELAY)
        run_singlezone(str(output_dir / dist.name), simtime,
                       overwrite=overwrite, RIa=dist, **STANDARD_PARAMS)
        plot_vice_onezone(str(output_dir / dist.name), fig=fig, axs=axs,
                          label=r'$\bar t={}$ Myr, $\sigma_t={}$ Myr'.format(
                              int(PEAKS[i]*1000), int(WIDTHS[i]*1000)),
                          style_kw={'linestyle': LINE_STYLE[i]},
                          marker_labels=(i==1),
        )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_prompt.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
