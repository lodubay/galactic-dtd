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
from track_and_mdf import setup_axes, plot_vice_onezone

# VICE one-zone model settings
DELAYS = [0.15, 0.04, 0.15, 0.04]
MODELS = ['exponential', 'exponential', 'powerlaw', 'powerlaw']
NRUNS = len(DELAYS)
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    tau_star=2.,
)

# Plot settings
LINE_STYLE = ['-', '--', '-', '--']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'delay_slope'

    fig, axs = setup_axes()

    for i in range(NRUNS):
        model = MODELS[i]
        delay = DELAYS[i]
        name = gen_name_from_params(model, delay)
        dtd_label = {
            'powerlaw': 'Power-Law',
            'exponential': 'Exponential'
        }
        label = rf'%s, $t_D={int(delay*1000)}$ Myr' % dtd_label[model]

        if overwrite:
            run(output_dir, i)
        else:
            try:
                history = vice.history(str(output_dir / name))
            except IOError:
                run(output_dir, i)

        if model == 'exponential':
            color = paultol.bright.colors[5]
            line_width = 1.5
            zorder = 1
        else:
            color = 'k'
            line_width = 1
            zorder = 10

        plot_vice_onezone(str(output_dir / name), fig=fig, axs=axs,
                          plot_kw={'label': label},
                          style_kw={
                              'color': color,
                              'linestyle': LINE_STYLE[i],
                              'linewidth': line_width,
                              'zorder': zorder},
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2, fontsize=7)
    fig.savefig(paths.figures / 'onezone_dtd_delay.pdf', dpi=300)
    plt.close()


def run(output_dir, i):
    """
    Set up and run the ith one-zone model.
    """
    sz = setup_single(output_dir, MODELS[i], DELAYS[i], **STANDARD_PARAMS)
    simtime = np.arange(0, END_TIME + DT, DT)
    sz.run(simtime, overwrite=True)


def setup_single(output_dir, model, delay, **kwargs):
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
    name = gen_name_from_params(model, delay)
    dtd = {
        'powerlaw': dtds.powerlaw(slope=-1.1, tmin=delay),
        'exponential': dtds.exponential(timescale=1.5, tmin=delay)
    }
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtd[model], delay=delay,
                         **kwargs)
    return sz


def gen_name_from_params(model, delay):
    """
    Generate singlezone output name from its distinguishing parameters.
    """
    name = f'{model}_delay{int(delay*1000):03d}'
    return name


if __name__ == '__main__':
    main()
