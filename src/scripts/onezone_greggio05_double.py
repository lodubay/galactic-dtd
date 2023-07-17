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

# DTD parameters
INTSTEP = 1e-4 # DTD integration timestep in Gyr
NSAMPLES = 200 # sets the integration precision; higher means longer runtime

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'greggio05'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    # Double-degenerate WIDE model
    dtd = dtds.greggio05_double('wide', dt=INTSTEP, nsamples=NSAMPLES)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Double Degenerate WIDE',
                      color=paultol.muted.colors[1],
                      marker_labels=True,
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1.5,
                           'zorder': 1},
                      logmdf=False
                      )
    
    # Double-degenerate CLOSE model
    dtd = dtds.greggio05_double('close', dt=INTSTEP, nsamples=NSAMPLES)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd, 
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label='Double Degenerate CLOSE',
                      color=paultol.muted.colors[3],
                      style_kw={
                           'linestyle': '-',
                           'linewidth': 1.5,
                           'zorder': 1},
                      logmdf=False
                      )
    
    dtd = dtds.plateau(width=1, slope=-1.1)
    sz = vice.singlezone(name=str(output_dir / dtd.name), 
                         RIa=dtd, 
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Plateau ($W=1$ Gyr, $\alpha=-1.1$)',
                      color=paultol.muted.colors[6],
                      style_kw={
                            'linestyle': '--',
                            'linewidth': 1,
                            'zorder': 10,
                      },
                      logmdf=False
                      )
    
    dtd = dtds.plateau(width=0.3, slope=-1.1)
    sz = vice.singlezone(name=str(output_dir / dtd.name), 
                          RIa=dtd, 
                          func=models.insideout(8, dt=dt), 
                          mode='sfr',
                          **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), fig=fig, axs=axs,
                      label=r'Plateau ($W=300$ Myr, $\alpha=-1.1$)',
                      color=paultol.muted.colors[7],
                      style_kw={
                            'linestyle': '--',
                            'linewidth': 1,
                            'zorder': 10,
                      },
                      logmdf=False)

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
    fig.savefig(paths.figures / 'onezone_greggio05_double.pdf', dpi=300)
    plt.close()    


if __name__ == '__main__':
    main()
