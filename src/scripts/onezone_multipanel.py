"""
This script generates plots for a multi-panel figure comparing the results
of several one-zone chemical evolution models with various DTD parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models, dtds
from _globals import END_TIME, ONEZONE_DEFAULTS
from colormaps import paultol
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)
    plot_powerlaw_slope()
    plot_exponential_timescale()
    plot_plateau_width()
    plot_minimum_delay()
    plot_all()


def plot_powerlaw_slope():
    slopes = [-0.8, -1.1,-1.4]
    line_styles = ['-', '--', ':']
    
    output_dir = paths.data / 'onezone' / 'powerlaw_slope'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    fig, axs = setup_axes(title='Power-law DTD')
    
    for i, slope in enumerate(slopes):
        dtd = dtds.powerlaw(slope=slope, tmin=ONEZONE_DEFAULTS['delay'])
        run_onezone(dtd, output_dir=output_dir)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs, 
                          linestyle=line_styles[i], 
                          color=styles.plaw['color'], 
                          label=rf'$\alpha={slope:.1f}$', 
                          marker_labels=(i==0))
    
    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_powerlaw_slope.pdf', dpi=300)
    plt.close()
    
    
def plot_exponential_timescale():
    timescales = [6, 3, 1.5]
    line_styles = ['-', '--', ':']
    
    output_dir = paths.data / 'onezone' / 'exponential_timescale'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(title='Exponential DTD')

    delay = ONEZONE_DEFAULTS['delay']
    for i, timescale in enumerate(timescales):
        dtd = dtds.exponential(timescale=timescale, tmin=delay)
        run_onezone(dtd, output_dir=output_dir)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs, 
                          linestyle=line_styles[i], 
                          color=styles.exp['color'], 
                          label=rf'$\tau={timescale:.1f}$ Gyr', 
                          marker_labels=(i==0))

    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_exponential_timescale.pdf', dpi=300)
    plt.close()


def plot_plateau_width():
    widths = [1., 0.3, 0.1]
    line_styles = ['-', '--', '-.']
    
    output_dir = paths.data / 'onezone' / 'plateau_width'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes(title='Plateau DTD')
    
    # Exponential for comparison
    delay = ONEZONE_DEFAULTS['delay']
    dtd = dtds.exponential(timescale=3, tmin=delay)
    run_onezone(dtd, output_dir=output_dir)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=fig, axs=axs,
                      linestyle=':', 
                      color=styles.plateau['color'], 
                      label=rf'Exponential ($\tau=3$ Gyr)', 
                      marker_labels=True, 
                      linewidth=1.5)

    slope = -1.1
    for i, width in enumerate(widths):
        dtd = dtds.plateau(width=width, slope=slope, tmin=delay)
        run_onezone(dtd, output_dir=output_dir)
        plot_vice_onezone(str(output_dir / dtd.name), 
                          fig=fig, axs=axs,
                          linestyle=line_styles[i], 
                          color=styles.exp['color'], 
                          label=r'$W={:1.1f}$ Gyr'.format(width))
    
    # Power law for comparison
    dtd = dtds.powerlaw(slope=slope, tmin=delay)
    run_onezone(dtd, output_dir=output_dir)
    plot_vice_onezone(str(output_dir / dtd.name), 
                      fig=fig, axs=axs, 
                      linestyle=':', 
                      color=styles.plaw['color'], 
                      label='No plateau')

    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_plateau_width.pdf', dpi=300)
    plt.close()


def plot_minimum_delay():
    delays = [0.15, 0.04]
    line_styles = ['--', '-']
    
    output_dir = paths.data / 'onezone' / 'minimum_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    for delay, ls in zip(delays, line_styles):
        ONEZONE_DEFAULTS['delay'] = delay
        distributions = [styles.exp, styles.plaw]
        for i, dtd in enumerate(distributions):
            if delay == delays[1]:
                label = dtd['label']
            else:
                label = None
            name = dtd['func'].name + '_delay{:03d}'.format(int(delay * 1000))
            run_onezone(dtd['func'], name=name, params={'delay': delay}, 
                        output_dir=output_dir)
            plot_vice_onezone(str(output_dir / name), 
                              fig=fig, axs=axs,
                              label=label, 
                              color=dtd['color'],
                              linestyle=ls,
                              zorder=i,
                              marker_labels=(i==0 and delay==delays[0])
                              )

    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    # Legend
    handles, labels = axs[0].get_legend_handles_labels()
    handles += [Line2D([], [], color='k', ls=ls, lw=1) for ls in line_styles]
    labels += [f'{int(1000*delay)} Myr minimum delay' for delay in delays]
    axs[0].legend(handles, labels, frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_minimum_delay.pdf', dpi=300)
    plt.close()


def plot_all():
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fig, axs = setup_axes()

    distributions = [styles.triple, styles.plateau_long, styles.exp, 
                     styles.plaw, styles.prompt]

    for i, dtd in enumerate(distributions):
        run_onezone(dtd['func'], output_dir=output_dir)
        plot_vice_onezone(str(output_dir / dtd['func'].name), 
                          fig=fig, axs=axs,
                          label=dtd['label'], 
                          color=dtd['color'],
                          marker_labels=(i==0),
                          )

    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)

    axs[0].legend(frameon=False, loc='lower left')
    fig.savefig(paths.figures / 'onezone_dtd.pdf', dpi=300)
    plt.close()


def run_onezone(dtd, params={}, output_dir='', name=''):
    """
    Run a one-zone model with the given DTD and parameters.
    
    Parameters
    ----------
    dtd : function
    params : dict, optional
        Other parameters passed to vice.singlezone. Overrides defaults in 
        ONEZONE_DEFAULTS.
    output_dir : str or Path, optional
    """
    for key in ONEZONE_DEFAULTS.keys():
        if key not in params.keys():
            params[key] = ONEZONE_DEFAULTS[key]
    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    if name == '':
        name = dtd.name
    sz = vice.singlezone(name=str(output_dir / name),
                         RIa=dtd,
                         func=models.insideout(8, dt=dt), 
                         mode='sfr',
                         **params)
    sz.run(simtime, overwrite=True)


if __name__ == '__main__':
    main()
