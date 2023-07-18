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
from _globals import END_TIME, ONEZONE_DEFAULTS, DT
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone

DELAYS = [0.15, 0.04]
LINE_STYLES = ['--', '-']

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'minimum_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    simtime = np.arange(0, END_TIME + DT, DT)

    labels = [r'Exponential ($\tau=3$ Gyr)',
              r'Plateau ($W=0.3$ Gyr)',
              r'Power-Law ($\alpha=-1.1$)',]
    colors = [paultol.vibrant.colors[i] for i in [0, 1, 4]]

    for delay, ls in zip(DELAYS, LINE_STYLES):
        ONEZONE_DEFAULTS['delay'] = delay
        distributions = [dtds.exponential(timescale=3, tmin=delay),
                         dtds.plateau(width=0.3, slope=-1.1, tmin=delay),
                         dtds.powerlaw(slope=-1.1, tmin=delay),]
        for i, dtd in enumerate(distributions):
            if delay == DELAYS[1]:
                label = labels[i]
            else:
                label = None
            name = dtd.name + '_delay{:03d}'.format(int(delay * 1000))
            sz = vice.singlezone(name=str(output_dir / name),
                                 RIa=dtd,
                                 func=models.insideout(8, dt=DT), 
                                 mode='sfr',
                                 **ONEZONE_DEFAULTS)
            sz.run(simtime, overwrite=True)
            plot_vice_onezone(str(output_dir / name), 
                              fig=fig, axs=axs,
                              label=label, color=colors[i],
                              style_kw={'linestyle': ls,
                                        'linewidth': 1,
                                        'zorder': 5-i+2*int(delay*10)},
                              marker_labels=(i==0 and delay==DELAYS[1])
                              )

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
