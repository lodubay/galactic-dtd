"""
This script plots abundance tracks from one-zone models with varying Type Ia
delay time distribution (DTD) and minimum delay time.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import vice
import paths
# from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_figure, plot_vice_onezone

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'yields'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_figure()

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    from multizone.src.yields import J21
    name = 'J21'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=1.5),
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), 
                      fig=fig, axs=axs,
                      label='Johnson+ (2021)', 
                      color='k',
                      linestyle='-',
                      # zorder=i,
                      # marker_labels=(i==0 and delay==DELAYS[0])
                      )

    from multizone.src.yields import C22
    name = 'C22'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=1.5),
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), 
                      fig=fig, axs=axs,
                      label='Conroy+ (2022)', 
                      color='r',
                      linestyle=':',
                      # zorder=i,
                      marker_labels=True
                      )

    from multizone.src.yields import W23
    name = 'W23'
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtds.exponential(timescale=1.5),
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / name), 
                      fig=fig, axs=axs,
                      label='Weinberg+ (2023)', 
                      color='b',
                      linestyle='--',
                      # zorder=i,
                      # marker_labels=(i==0 and delay==DELAYS[0])
                      )


    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    # Legend
    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    fig.savefig(paths.figures / 'onezone_yields.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
