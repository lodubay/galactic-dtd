"""
This script plots abundance tracks from one-zone models with a broken power-law
DTD, varying the length of the initial "plateau" in the Ia rate.
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

PLATEAUS = [1., 0.3, 0.1] # Gyr
SLOPE = -1.1
# Plot settings
LINE_STYLE = ['-', '--', '-.']
LOG_MDF = False

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'plateau_width'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    dt = ONEZONE_DEFAULTS['dt']
    delay = ONEZONE_DEFAULTS['delay']
    simtime = np.arange(0, END_TIME + dt, dt)

    fig, axs = setup_axes(logmdf=LOG_MDF, title='Plateau DTD')

    # Plot exponentials for reference
    # for tau, ls in zip([1.5, 3], ['--', '-']):
    tau = 3
    exp = dtds.exponential(timescale=tau, tmin=delay)
    sz = vice.singlezone(name=str(output_dir / exp.name),
                         RIa=exp,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    
    plot_vice_onezone(str(output_dir / exp.name),
                      fig=fig, axs=axs,
                      label=rf'Exponential ($\tau={tau:.01f}$ Gyr)',
                      color=paultol.bright.colors[0],
                      style_kw={'linestyle': ':',
                                'linewidth': 1.5,
                                'zorder': 1},
                      logmdf=LOG_MDF
                      )

    for i, plateau in enumerate(PLATEAUS):
        if plateau >= 1:
            label = f'{plateau:.01f} Gyr plateau'
        else:
            label = f'{int(plateau*1000)} Myr plateau'

        dtd = dtds.plateau(width=plateau, slope=SLOPE, tmin=delay)
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)

        plot_vice_onezone(str(output_dir / dtd.name),
                          fig=fig, axs=axs,
                          label=r'$W={:1.1f}$ Gyr'.format(plateau), 
                          color=paultol.bright.colors[1],
                          style_kw={
                              'linestyle': LINE_STYLE[i],
                              'linewidth': 1,
                              'zorder': 9},
                          marker_labels=(i==0),
                          logmdf=LOG_MDF
                          )

    # Plot standard power-law for reference
    plaw = dtds.powerlaw(slope=SLOPE, tmin=delay)
    sz = vice.singlezone(name=str(output_dir / plaw.name),
                         RIa=plaw,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / plaw.name),
                      fig=fig, axs=axs,
                      label='No plateau', color='k',
                      style_kw={'linestyle': ':',
                                'linewidth': 1,
                                'zorder': 2},
                      logmdf=LOG_MDF
                      )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.4))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_plateau_width.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()