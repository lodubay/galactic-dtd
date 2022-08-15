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
from utils import run_singlezone

# VICE one-zone model settings
PLATEAUS = [0.1, 0.3, 1.] # Gyr
DELAY = 0.04 # minimum delay time in Gyr
DT = 0.01
SLOPE = -1.1
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
LINE_STYLE = [':', '-.', '-']
LOG_MDF = False

def main(overwrite=False):
    output_dir = paths.data / 'onezone' / 'plateau'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    simtime = np.arange(0, END_TIME + DT, DT)

    fig, axs = setup_axes(logmdf=LOG_MDF)

    # Plot standard power-law for reference
    dist = dtds.powerlaw(slope=SLOPE, tmin=DELAY)
    run_singlezone(str(output_dir / dist.name), simtime, overwrite=overwrite,
                   RIa=dist, **STANDARD_PARAMS)
    plot_vice_onezone(str(output_dir / dist.name),
                      fig=fig, axs=axs,
                      label='No plateau', color='k',
                      style_kw={'linestyle': '-',
                                'linewidth': 1,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )

    for i, plateau in enumerate(PLATEAUS):
        if plateau >= 1:
            label = f'{plateau:.01f} Gyr plateau'
        else:
            label = f'{int(plateau*1000)} Myr plateau'

        dist = dtds.plateau(width=plateau, slope=SLOPE, tmin=DELAY)
        run_singlezone(str(output_dir / dist.name), simtime, overwrite=overwrite,
                       RIa=dist, **STANDARD_PARAMS)

        plot_vice_onezone(str(output_dir / dist.name),
                          fig=fig, axs=axs,
                          label=label, color=paultol.bright.colors[1],
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1,
                              'zorder': 10},
                          logmdf=LOG_MDF
                          )

    # Plot exponentials for reference
    for tau, ls in zip([1.5, 3], ['--', '-']):
        dist = dtds.exponential(timescale=tau, tmin=DELAY)
        run_singlezone(str(output_dir / dist.name), simtime, overwrite=overwrite,
                       RIa=dist, **STANDARD_PARAMS)

        plot_vice_onezone(str(output_dir / dist.name),
                          fig=fig, axs=axs,
                          label=rf'Exponential ($\tau={tau:.01f}$ Gyr)',
                          color=paultol.bright.colors[0],
                          style_kw={'linestyle': ls,
                                    'linewidth': 1.5,
                                    'zorder': 1},
                          marker_labels=(tau==3),
                          logmdf=LOG_MDF
                          )

    # Adjust axis limits
    axs[0].set_xlim((-2.5, 0.3))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.5, fontsize=7)
    fig.savefig(paths.figures / 'onezone_plateau.pdf', dpi=300)
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
