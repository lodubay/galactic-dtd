"""
This script plots abundance tracks for one-zone models with varying minimum
Type Ia delay times and star formation efficiency timescales.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import models, dtds
from migration.src._globals import END_TIME
from colormaps import paultol
from tri_panel_ism_tracks import setup_axes

# Settings
MINIMUM_DELAY = [0.10, 0.05, 0.20, 0.10, 0.10] # Gyr
TAU_STAR = [2.0, 2.0, 2.0, 1.0, 4.0] # Gyr
LINE_STYLE = ['-', ':', '--', '-', '-']
COLOR = ['k', 'k', 'k', paultol.highcontrast.colors[2],
         paultol.highcontrast.colors[1]]
NRUNS = len(MINIMUM_DELAY)
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    RIa=dtds.powerlaw(),
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
)

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'delay_taustar'

    fig, axs = setup_axes(felim=(-2., 0.2), tight_layout=True)

    for i in range(NRUNS):
        delay = MINIMUM_DELAY[i]
        tau_star = TAU_STAR[i]
        name = gen_name_from_params(delay=delay, tau_star=tau_star)
        label = rf'$t_D={int(delay*1000)}$ Myr, $\tau_*={int(tau_star)}$ Gyr'

        if overwrite:
            run(output_dir, i)
        else:
            try:
                history = vice.history(str(output_dir / name))
            except IOError:
                run(output_dir, i)
        history = vice.history(str(output_dir / name))

        if tau_star != 2.0:
            line_width = 2
            zorder = 1
        else:
            line_width = 1.5
            zorder = 10

        axs[0].plot(history['time'], history['[fe/h]'], label=label,
                    color=COLOR[i], ls=LINE_STYLE[i], lw=line_width,
                    zorder=zorder)
        axs[1].plot(history['time'], history['[o/fe]'], color=COLOR[i],
                    ls=LINE_STYLE[i], lw=line_width, zorder=zorder)
        axs[2].plot(history['[fe/h]'], history['[o/fe]'],
                    color=COLOR[i], ls=LINE_STYLE[i], lw=line_width,
                    zorder=zorder)

    axs[0].legend(frameon=False)
    fig.savefig(paths.figures / 'onezone_delay_taustar.png', dpi=300)
    plt.close()


def run(output_dir, i):
    """
    Set up and run the ith one-zone model.
    """
    sz = setup_single(output_dir,
                      delay=MINIMUM_DELAY[i], tau_star=TAU_STAR[i],
                      **STANDARD_PARAMS)
    simtime = np.arange(0, END_TIME + DT, DT)
    sz.run(simtime, overwrite=True)


def setup_single(output_dir, delay=0.1, tau_star=2., **kwargs):
    """
    Setup a one-zone model with a given minimum Ia delay time and SFE timescale.

    Parameters
    ----------
    output_dir : Path
        Parent directory for VICE outputs
    delay : float, optional
        Minimum Type Ia delay time in Gyr. The default is 0.1 Gyr.
    tau_star : float, optional
        Star formation efficiency timescale in Gyr. The default is 2 Gyr.
    Other keyword arguments are passed to vice.singlezone

    Returns
    -------
    sz : vice.singlezone object
    """
    name = gen_name_from_params(delay=delay, tau_star=tau_star)
    sz = vice.singlezone(name=str(output_dir / name),
                         delay=delay, tau_star=tau_star,
                         **kwargs)
    return sz


def gen_name_from_params(delay=0.1, tau_star=2.0):
    """
    Generate singlezone output name from its distinguishing parameters.
    """
    name = 'delay{:03d}_taustar{:02d}'.format(int(delay*1000), int(tau_star*10))
    return name



if __name__ == '__main__':
    main()
