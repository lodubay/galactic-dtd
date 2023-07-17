"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
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

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    delay = ONEZONE_DEFAULTS['delay']
    distributions = [dtds.triple(tmin=delay),
                     dtds.plateau(width=1., slope=-1.1, tmin=delay),
                     dtds.exponential(timescale=1.5, tmin=delay),
                     dtds.powerlaw(slope=-1.1, tmin=delay),
                     dtds.prompt(peak=0.05, stdev=0.015, timescale=3, tmin=delay),
                     ]
    labels = [r'Triple-system',
              r'Plateau ($W=1$ Gyr)',
              r'Exponential ($\tau=1.5$ Gyr)',
              r'Power-law ($\alpha=-1.1$)',
              r'Prompt',]
    colors = [paultol.vibrant.colors[i] for i in [5, 0, 1, 4, 2]]
    line_styles = ['-', '-.', '--', ':', ':']

    for i, dtd in enumerate(distributions):
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          label=labels[i], 
                          color=colors[i],
                          style_kw={'linestyle': '-',
                                    'linewidth': 1},
                          marker_labels=(i==0),
                          # mdf_style='curve'
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.4))
    axs[0].set_ylim((-0.1, 0.52))
    axs[1].set_ylim((0, None))
    axs[2].set_xlim((0, None))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
