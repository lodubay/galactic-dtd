"""
This script plots abundance tracks from one-zone models with a power-law DTD
with varying slope.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
from multizone.src import models, dtds
from _globals import END_TIME
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from utils import run_singlezone

# VICE one-zone model settings
SLOPES = [-0.8, -1.1,-1.4]
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
COLOR = 'k'
LOG_MDF = False

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'powerlaw_slope'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=LOG_MDF, title='Power-law DTD')

    simtime = np.arange(0, END_TIME + DT, DT)

    for i, slope in enumerate(SLOPES):
        dtd = dtds.powerlaw(slope=slope, tmin=DELAY)
        run_singlezone(str(output_dir / dtd.name), simtime,
                       overwrite=True, RIa=dtd, **STANDARD_PARAMS)

        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          label=rf'$\alpha={slope:.1f}$',
                          color=COLOR,
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1},
                          logmdf=LOG_MDF,
                          marker_labels=(i==0)
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.0, 0.4))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_powerlaw_slope.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
