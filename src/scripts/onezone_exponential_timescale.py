"""
This script plots abundance tracks from one-zone models with an exponential DTD
with varying timescale.
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

TIMESCALES = [6, 3, 1.5]
LINE_STYLE = ['-', '--', ':']
COLOR = paultol.bright.colors[0]
LOG_MDF = False

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'exponential_timescale'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=LOG_MDF, title='Exponential DTD')

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    for i, timescale in enumerate(TIMESCALES):
        dtd = dtds.exponential(timescale=timescale, 
                               tmin=ONEZONE_DEFAULTS['delay'])
        # Run one-zone model
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        # Plot one-zone models
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          label=rf'$\tau={timescale:.1f}$ Gyr',
                          color=paultol.bright.colors[0],
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1},
                          marker_labels=(i==0),
                          logmdf=LOG_MDF
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.4))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_exponential_timescale.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
