"""
This script plots abundance tracks from one-zone models with a broken power-law
DTD, varying the length of the initial "plateau" in the Ia rate.
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
PLATEAUS = [0.1, 0.2, 0.5, 1.] # Myr
DELAY = 0.04 # minimum delay time in Gyr
DT = 0.01
STANDARD_PARAMS = dict(
    func=models.insideout(8, dt=DT),
    mode='sfr',
    elements=('fe', 'o'),
    dt=DT,
    recycling='continuous',
    eta=2.5,
    tau_star=2.,
    delay=DELAY,
)

# Plot settings
LINE_STYLE = [':', '-.', '--', '-']

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'powerlaw_plateau'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    simtime = np.arange(0, END_TIME + DT, DT)

    fig, axs = setup_axes()

    # Plot standard power-law for reference
    sz = vice.singlezone(name=str(output_dir / 'powerlaw'),
                         RIa=dtds.powerlaw(slope=-1.1),
                         **STANDARD_PARAMS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / 'powerlaw'),
                      fig=fig, axs=axs,
                      plot_kw={'label': 'No plateau'},
                      style_kw={'color': 'k',
                                'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 1},
                      )

    for i, plateau in enumerate(PLATEAUS):
        name = gen_name_from_params(plateau)
        if plateau >= 1:
            label = f'{plateau:.01f} Gyr plateau'
        else:
            label = f'{int(plateau*1000)} Myr plateau'

        if overwrite:
            sz = vice.singlezone(name=str(output_dir / name),
                                 RIa=dtds.powerlaw_broken(tsplit=plateau+DELAY),
                                 **STANDARD_PARAMS)
            sz.run(simtime, overwrite=True)
        else:
            try:
                history = vice.history(str(output_dir / name))
            except IOError:
                sz = vice.singlezone(name=str(output_dir / name),
                                     RIa=dtds.powerlaw_broken(tsplit=plateau+DELAY),
                                     **STANDARD_PARAMS)
                sz.run(simtime, overwrite=True)

        plot_vice_onezone(str(output_dir / name),
                          fig=fig, axs=axs,
                          plot_kw={'label': label},
                          style_kw={
                              'color': paultol.bright.colors[1],
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1,
                              'zorder': 10},
                          )

    # Plot exponentials for reference
    for tau, ls in zip([1.5, 3], ['--', '-']):
        name = f'exponential{int(tau*10):02d}'
        sz = vice.singlezone(name=str(output_dir / name),
                             RIa=dtds.exponential(timescale=tau),
                             **STANDARD_PARAMS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / name),
                          fig=fig, axs=axs,
                          plot_kw={'label': rf'Exponential ($\tau={tau:.01f}$ Gyr)'},
                          style_kw={'color': paultol.bright.colors[0],
                                    'linestyle': ls,
                                    'linewidth': 1.5,
                                    'zorder': 1},
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.2))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.5, fontsize=7)
    fig.savefig(paths.figures / 'onezone_powerlaw_plateau.pdf', dpi=300)
    plt.close()


# def run(output_dir, i):
#     """
#     Set up and run the ith one-zone model.
#     """
#     sz = setup_single(output_dir, MODELS[i], DELAYS[i], **STANDARD_PARAMS)
#     simtime = np.arange(0, END_TIME + DT, DT)
#     sz.run(simtime, overwrite=True)


# def setup_single(output_dir, model, delay, **kwargs):
#     """
#     Setup a one-zone model with a given minimum Ia delay time and SFE timescale.

#     Parameters
#     ----------
#     output_dir : Path
#         Parent directory for VICE outputs
#     delay : float, optional
#         Minimum Type Ia delay time in Gyr. The default is 0.1 Gyr.
#     tau_star : float, optional
#         Star formation efficiency timescale in Gyr. The default is 2 Gyr.
#     Other keyword arguments are passed to vice.singlezone

#     Returns
#     -------
#     sz : vice.singlezone object
#     """
#     name = gen_name_from_params(model, delay)
#     dtd = {
#         'powerlaw': dtds.powerlaw(slope=-1.1, tmin=delay),
#         'exponential': dtds.exponential(timescale=1.5, tmin=delay)
#     }
#     sz = vice.singlezone(name=str(output_dir / name),
#                          RIa=dtd[model], delay=delay,
#                          **kwargs)
#     return sz


def gen_name_from_params(plateau):
    """
    Generate singlezone output name from its distinguishing parameters.
    """
    name = f'plateau{int(plateau*1000):03d}'
    return name


if __name__ == '__main__':
    main()
