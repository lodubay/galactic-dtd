"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_figure, plot_vice_onezone
from delay_time_distributions import styles

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_figure()

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    delay = ONEZONE_DEFAULTS['delay']
    distributions = [styles.gauss_narrow, styles.gauss_wide, styles.triple, 
                     styles.plateau_long, styles.exp, 
                     styles.plaw, styles.prompt]

    for i, dtd in enumerate(distributions):
        sz = vice.singlezone(name=str(output_dir / dtd['func'].name),
                             RIa=dtd['func'],
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd['func'].name), 
                          fig=fig, axs=axs,
                          label=dtd['label'], 
                          color=dtd['color'],
                          linestyle=dtd['line'],
                          marker_labels=(i==0),
                          )

    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.8)
    fig.savefig(paths.extra / 'onezone_strolger05', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='onezone_dtd.py',
        description='Plot abundance tracks from one-zone models with varying DTD')
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
