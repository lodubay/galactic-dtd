"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD) and minimum delay time.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

DELAYS = [0.15, 0.04]
LINE_STYLES = ['--', '-']

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'minimum_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    for delay, ls in zip(DELAYS, LINE_STYLES):
        ONEZONE_DEFAULTS['delay'] = delay
        distributions = [styles.exp, styles.plaw]
        for i, dtd in enumerate(distributions):
            if delay == DELAYS[1]:
                label = dtd['label']
            else:
                label = None
            name = dtd['func'].name + '_delay{:03d}'.format(int(delay * 1000))
            sz = vice.singlezone(name=str(output_dir / name),
                                 RIa=dtd['func'],
                                 func=models.insideout(8, dt=dt), 
                                 mode='sfr',
                                 **ONEZONE_DEFAULTS)
            sz.run(simtime, overwrite=True)
            plot_vice_onezone(str(output_dir / name), 
                              fig=fig, axs=axs,
                              label=label, 
                              color=dtd['color'],
                              linestyle=ls,
                              zorder=i,
                              marker_labels=(i==0 and delay==DELAYS[0])
                              )


    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    # Legend
    handles, labels = axs[0].get_legend_handles_labels()
    handles += [Line2D([], [], color='k', ls=ls, lw=1) for ls in LINE_STYLES]
    labels += [f'{int(1000*delay)} Myr minimum delay' for delay in DELAYS]
    axs[0].legend(handles, labels, frameon=False, loc='lower left',
                  handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_minimum_delay.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
