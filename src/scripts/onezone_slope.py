"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD).
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
SLOPE = [-0.8, -1.1, -1.4]
LINE_STYLE = [':', '-', '--']
COLOR = ['k', 'k', 'k']
NRUNS = len(SLOPE)
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    delay=0.04,
    tau_star=2.,
)

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'slope'

    fig, axs = setup_axes(felim=(-2., 0.2), tight_layout=True)

    for i in range(NRUNS):
        slope = SLOPE[i]
        name = gen_name_from_params(slope=slope)
        label = rf'$\alpha={slope:.2f}$'

        if overwrite:
            run(output_dir, i)
        else:
            try:
                history = vice.history(str(output_dir / name))
            except IOError:
                run(output_dir, i)
        history = vice.history(str(output_dir / name))

        axs[0].plot(history['time'], history['[fe/h]'], label=label,
                    color=COLOR[i], ls=LINE_STYLE[i])
        axs[1].plot(history['time'], history['[o/fe]'],
                    color=COLOR[i], ls=LINE_STYLE[i])
        axs[2].plot(history['[fe/h]'], history['[o/fe]'],
                    color=COLOR[i], ls=LINE_STYLE[i])

    axs[0].legend(frameon=False)
    fig.savefig(paths.figures / 'onezone_slope.png', dpi=300)
    plt.close()


def run(output_dir, i):
    """
    Set up and run the ith one-zone model.
    """
    sz = setup_single(output_dir, slope=SLOPE[i], **STANDARD_PARAMS)
    simtime = np.arange(0, END_TIME + DT, DT)
    sz.run(simtime, overwrite=True)


def setup_single(output_dir, slope=-1.1, **kwargs):
    """
    Setup a one-zone model with a given minimum Ia delay time and SFE timescale.

    Parameters
    ----------
    output_dir : Path
        Parent directory for VICE outputs
    slope : float, optional
        Slope of the power-law delay time distribution. The default is -1.1.
    Other keyword arguments are passed to vice.singlezone

    Returns
    -------
    sz : vice.singlezone object
    """
    name = gen_name_from_params(slope=slope)
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.powerlaw(slope=slope),
                         **kwargs)
    return sz


def gen_name_from_params(slope=-1.1):
    """
    Generate singlezone output name from its distinguishing parameters.
    """
    name = 'slope{:03d}'.format(int(-10*slope))
    return name



if __name__ == '__main__':
    main()

