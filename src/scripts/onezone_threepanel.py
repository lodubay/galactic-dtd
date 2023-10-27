"""
Make a three-panel plot of onezone chemical evolution tracks for the power-law,
exponential, and plateau DTD models illustrating the effect of changing the
key parameter for each model.
"""

import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS, TWO_COLUMN_WIDTH
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig = plt.figure(figsize=(TWO_COLUMN_WIDTH, 0.36*TWO_COLUMN_WIDTH))
    gs = fig.add_gridspec(7, 22, wspace=0.)
    # subfigs = fig.subfigures(1, 3, wspace=0.)
    subfigs = [fig.add_subfigure(gs[:,i:i+w]) for i, w in zip((0, 8, 15), (8, 7, 7))]
    plaw_axs = powerlaw_slope(subfigs[0])
    exp_axs = exponential_timescale(subfigs[1])
    plat_axs = plateau_width(subfigs[2])
    plt.subplots_adjust(bottom=0.13, top=0.98, left=0.16, right=0.98, wspace=0.5)
    fig.savefig(paths.figures / 'onezone_threepanel.pdf', dpi=300)
    plt.close()


def powerlaw_slope(subfig, slopes=[-0.8, -1.1, -1.4], 
                   line_styles=['-', '--', ':']):
    output_dir = paths.data / 'onezone' / 'powerlaw_slope'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    axs = setup_axes(subfig, title='Power-law DTD', xlim=(-2.1, 0.6))
    # subfig.gridspec.update(right=1.)

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    
    # Exponential DTD for reference
    dtd = dtds.exponential(timescale=3, tmin=ONEZONE_DEFAULTS['delay'])
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=subfig, axs=axs, 
                      linestyle='-', 
                      linewidth=2,
                      color='#bbbbbb', 
                      label='Exp.\n($\\tau=3$ Gyr)', 
                      marker_labels=True,
                      zorder=1)

    for i, slope in enumerate(slopes):
        dtd = dtds.powerlaw(slope=slope, tmin=ONEZONE_DEFAULTS['delay'])
        # Run one-zone model
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=subfig, axs=axs, 
                          linestyle=line_styles[i], 
                          color=styles.plaw['color'], 
                          label=rf'$\alpha={slope:.1f}$', 
                          marker_labels=False)

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')

    return axs


def exponential_timescale(subfig, timescales=[6, 3, 1.5], 
                          line_styles=['-', '--', ':']):
    output_dir = paths.data / 'onezone' / 'exponential_timescale'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    axs = setup_axes(subfig, title='Exponential DTD', ylabel=False, xlim=(-2.1, 0.6))

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)

    for i, timescale in enumerate(timescales):
        dtd = dtds.exponential(timescale=timescale, 
                               tmin=ONEZONE_DEFAULTS['delay'])
        # Run one-zone model
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        # Plot one-zone models
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=subfig, axs=axs, 
                          linestyle=line_styles[i], 
                          color=styles.exp['color'], 
                          label=rf'$\tau={timescale:.1f}$ Gyr', 
                          marker_labels=(i==0))

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    
    return axs


def plateau_width(subfig, plateau_widths=[1., 0.3, 0.1], 
                  line_styles=['-', '--', '-.'], slope=-1.1):
    output_dir = paths.data / 'onezone' / 'plateau_width'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    dt = ONEZONE_DEFAULTS['dt']
    delay = ONEZONE_DEFAULTS['delay']
    simtime = np.arange(0, END_TIME + dt, dt)

    axs = setup_axes(subfig, title='Plateau DTD', ylabel=False, xlim=(-2.1, 0.6))

    # Plot exponentials for reference
    dtd = dtds.exponential(timescale=3, tmin=delay)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=subfig, axs=axs,
                      linestyle='-', 
                      linewidth=2,
                      color='#bbbbbb', 
                      label='Exp.\n($\\tau=3$ Gyr)', 
                      marker_labels=True, 
                      zorder=1)

    for i, width in enumerate(plateau_widths):
        dtd = dtds.plateau(width=width, slope=slope, tmin=delay)
        sz = vice.singlezone(name=str(output_dir / dtd.name),
                             RIa=dtd,
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=subfig, axs=axs,
                          linestyle=line_styles[i], 
                          color=styles.plateau['color'], 
                          label=r'$W={:1.1f}$ Gyr'.format(width))

    # Plot standard power-law for reference
    dtd = dtds.powerlaw(slope=slope, tmin=delay)
    sz = vice.singlezone(name=str(output_dir / dtd.name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=subfig, axs=axs, 
                      linestyle=':', 
                      color=styles.plaw['color'], 
                      label='No plateau')

    # Adjust axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    
    return axs


if __name__ == '__main__':
    main()
