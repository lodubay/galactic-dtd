import sys
import glob
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

# Settings
MINIMUM_DELAY = [0.05, 0.10, 0.20, 0.10, 0.10] # Gyr
TAU_STAR = [2.0, 2.0, 2.0, 1.0, 4.0] # Gyr
LINE_STYLE = [':', '-', '--', '-', '-']
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

    fig, ax = plt.subplots()

    for i in range(NRUNS):
        delay = MINIMUM_DELAY[i]
        tau_star = TAU_STAR[i]
        name = gen_name_from_params(delay=delay, tau_star=tau_star)
        label = rf'$t_D={delay}$ Myr, $\tau_*={tau_star}$ Gyr'

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

        ax.plot(history['[fe/h]'], history['[o/fe]'], label=label,
                color=COLOR[i], ls=LINE_STYLE[i], lw=line_width, zorder=zorder)

    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('[O/Fe]')
    ax.legend()
    fig.savefig(paths.figures / 'onezone_delay_taustar.pdf')
    plt.close()


def run(output_dir, i):
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
