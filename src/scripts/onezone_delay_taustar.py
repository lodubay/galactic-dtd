"""
This script plots abundance tracks for one-zone models with varying minimum
Type Ia delay times and star formation efficiency timescales.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, DT, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone

# One-zone model settings
MINIMUM_DELAY = [0.16, 0.08, 0.04, 0.08, 0.08] # Gyr
TAU_STAR = [2.0, 2.0, 2.0, 1.0, 4.0] # Gyr
NRUNS = len(MINIMUM_DELAY)

# Plot settings
LINE_STYLE = ['--', '-', ':', '-', '-']
COLOR = ['k', 'k', 'k', paultol.highcontrast.colors[2],
         paultol.highcontrast.colors[1]]

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'delay_taustar'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    simtime = np.arange(0, END_TIME + DT, DT)

    fig, axs = setup_axes(logmdf=False)

    for i in range(NRUNS):
        delay = MINIMUM_DELAY[i]
        tau_star = TAU_STAR[i]
        name = gen_name_from_params(delay=delay, tau_star=tau_star)
        label = rf'$t_D={int(delay*1000)}$ Myr, $\tau_*={int(tau_star)}$ Gyr'

        ONEZONE_DEFAULTS['tau_star'] = tau_star
        ONEZONE_DEFAULTS['delay'] = delay
        sz = vice.singlezone(name=str(output_dir / name),
                             RIa=dtds.powerlaw(slope=-1.1, tmin=delay),
                             func=models.insideout(8, dt=DT), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)

        if tau_star != 2.0:
            line_width = 1.5
            zorder = 1
        else:
            line_width = 1
            zorder = 10

        if delay == 0.16 and tau_star == 2:
            marker_labels = True
        else:
            marker_labels = False

        plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                          label=label, color=COLOR[i],
                          marker_labels=marker_labels,
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': line_width,
                              'zorder': zorder},
                          logmdf=False
                          )

    # Adjust axis limits
    # axs[0].set_xlim((-2.7, 0.3))
    # axs[0].set_ylim((-0.1, 0.54))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_delay_taustar.pdf', dpi=300)
    plt.close()


def gen_name_from_params(delay=0.1, tau_star=2.0):
    """
    Generate singlezone output name from its distinguishing parameters.
    """
    name = 'delay{:03d}_taustar{:02d}'.format(int(delay*1000), int(tau_star*10))
    return name


if __name__ == '__main__':
    main()
