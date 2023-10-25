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
from _globals import END_TIME, ONEZONE_DEFAULTS, TWO_COLUMN_WIDTH
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

SLOPES = [-0.8, -1.1,-1.4]
LINE_STYLES = ['-', '--', ':']

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'powerlaw_slope'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(width=0.4*TWO_COLUMN_WIDTH, title='Power-law DTD')

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    
    # Exponential DTD for reference
    dtd = dtds.exponential(timescale=3, tmin=ONEZONE_DEFAULTS['delay'])
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=fig, axs=axs, 
                      linestyle='-', 
                      linewidth=2,
                      color='#bbbbbb', 
                      label='Exponential\n($\\tau=3$ Gyr)', 
                      marker_labels=True,
                      zorder=1)

    for i, slope in enumerate(SLOPES):
        dtd = dtds.powerlaw(slope=slope, tmin=ONEZONE_DEFAULTS['delay'])
        # Run one-zone model
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs, 
                          linestyle=LINE_STYLES[i], 
                          color=styles.plaw['color'], 
                          label=rf'$\alpha={slope:.1f}$', 
                          marker_labels=False)

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_powerlaw_slope.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
