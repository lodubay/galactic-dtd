"""
This script plots abundance tracks from one-zone models with a broken power-law
DTD, varying the length of the initial "plateau" in the Ia rate.
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

PLATEAU_WIDTHS = [1., 0.3, 0.1] # Gyr
SLOPE = -1.1
# Plot settings
LINE_STYLES = ['-', '--', '-.']

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'plateau_width'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    dt = ONEZONE_DEFAULTS['dt']
    delay = ONEZONE_DEFAULTS['delay']
    simtime = np.arange(0, END_TIME + dt, dt)

    fig, axs = setup_axes(width=0.4*TWO_COLUMN_WIDTH, title='Plateau DTD')

    # Plot exponentials for reference
    dtd = dtds.exponential(timescale=3, tmin=delay)
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

    for i, width in enumerate(PLATEAU_WIDTHS):
        dtd = dtds.plateau(width=width, slope=SLOPE, tmin=delay)
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          linestyle=LINE_STYLES[i], 
                          color=styles.plateau['color'], 
                          label=r'$W={:1.1f}$ Gyr'.format(width))

    # Plot standard power-law for reference
    dtd = dtds.powerlaw(slope=SLOPE, tmin=delay)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=fig, axs=axs, 
                      linestyle=':', 
                      color=styles.plaw['color'], 
                      label='No plateau')

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)

    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_plateau_width.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
