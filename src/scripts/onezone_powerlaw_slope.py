"""
This script plots abundance tracks from one-zone models with a power-law DTD
with varying slope.
"""

import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone

SLOPES = [-0.8, -1.1,-1.4]
LINE_STYLES = ['-', '--', ':']
COLOR = 'k'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'powerlaw_slope'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(title='Power-law DTD')

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    for i, slope in enumerate(SLOPES):
        dtd = dtds.powerlaw(slope=slope, tmin=ONEZONE_DEFAULTS['delay'])
        # Run one-zone model
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        # Plot one-zone model
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          label=rf'$\alpha={slope:.1f}$',
                          color=COLOR,
                          linestyle=LINE_STYLES[i],
                          marker_labels=(i==0)
                          )

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)

    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_powerlaw_slope.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
