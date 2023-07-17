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
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'greggio05'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(logmdf=False)

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    # Single-degenerate
    dtd = dtds.greggio05_single()
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, 
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Single Degenerate',
                      color=paultol.muted.colors[5],
                      style_kw={
                          'linestyle': '-',
                          'linewidth': 1.5,
                          'zorder': 1},
                      logmdf=False
                      )

    dtd = dtds.exponential(timescale=1.5)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, 
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Exponential ($\tau=1.5$ Gyr)',
                      color=paultol.muted.colors[0],
                      style_kw={
                            'linestyle': '--',
                            'linewidth': 1,
                            'zorder': 10,
                      }, 
                      logmdf=False,
                      marker_labels=True
                      )

    dtd = dtds.powerlaw(slope=-1.1)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, 
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Power-Law ($\alpha=-1.1$)',
                      color='k',
                      style_kw={
                           'linestyle': ':',
                           'linewidth': 1,
                           'zorder': 1,
                      }, 
                      logmdf=False
                      )

    # Adjust axis limits
    axs[0].set_xlim((-2.1, 0.4))
    axs[0].set_ylim((-0.1, 0.52))

    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_greggio05_single.pdf', dpi=300)
    plt.close()    


if __name__ == '__main__':
    main()
