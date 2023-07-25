"""
This script generates the two-panel figure about the Greggio (2005) analytical
DTDs. The first plot shows the DTDs as a function of time, and the second
runs one-zone models and plots the outputs in abundance space. Reference this 
script in the `\script{}` command when plotting both figures side-by-side in 
the manuscript.
"""
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src import models, dtds
from _globals import END_TIME, MIN_RIA_DELAY, ONEZONE_DEFAULTS
from colormaps import paultol

INTSTEP = 1e-4 # plotting timestep in Gyr
NSAMPLES = 200 # sets the integration precision; higher means longer runtime

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    plot_dtd()
    plot_onezone()


def plot_dtd():
    from delay_time_distributions import setup_axes
    fig, ax = setup_axes()
    tarr = np.logspace(np.log10(MIN_RIA_DELAY), np.log10(END_TIME), 
                       num=NSAMPLES)
    
    for dtd in styles.distlist:
        func = dtd['func']
        ax.plot(tarr, dtd['offset'] * 1e9 * np.array([func(t) for t in tarr]), 
                label=dtd['label'], c=dtd['color'], ls=dtd['linestyle'], lw=1)
    
    ax.set_ylim((3e-3, 10))
    ax.legend(loc='lower left', frameon=False, bbox_to_anchor=(0.1, 0.01))
    ax.set_ylabel(r'Normalized SN Ia rate $\times$ factor')
    
    fig.savefig(paths.figures / 'dtd_analytical.pdf', dpi=300)
    plt.close()


def plot_onezone():
    from track_and_mdf import setup_axes, plot_vice_onezone
    from multizone.src.yields import J21
    
    output_dir = paths.data / 'onezone' / 'analytical_dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    fig, axs = setup_axes(xlim=(-1.9, 0.6))
    
    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    params = ONEZONE_DEFAULTS
    
    for i, dtd in enumerate(styles.distlist):
        params['eta'] = dtd['offset']
        sz = vice.singlezone(name=str(output_dir / dtd['func'].name),
                             RIa=dtd['func'],
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **params)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd['func'].name), 
                          fig=fig, axs=axs,
                          color=dtd['color'],
                          marker_labels=(i==4),
                          linestyle=dtd['linestyle']
                          )
    
    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    fig.savefig(paths.figures / 'onezone_analytical_dtd.pdf', dpi=300)
    plt.close()  


class styles:
    """Plot styling for different DTD models."""
    # Greggio (2005) DTDs
    sd = {
        'func': dtds.greggio05_single(),
        'label': r'Single Degenerate $\times6$',
        'color': paultol.muted.colors[5],
        'linestyle': '-',
        'offset': 4,
    }
    dd_wide = {
        'func': dtds.greggio05_double('wide', dt=INTSTEP, nsamples=NSAMPLES,
                                      progress=True),
        'label': 'Double Degenerate (WIDE)',
        'color': paultol.muted.colors[1],
        'linestyle': '-',
        'offset': 1,
        'eta': 1,
    }
    dd_close = {
        'func': dtds.greggio05_double('close', dt=INTSTEP, nsamples=NSAMPLES,
                                      progress=True),
        'label': r'Double Degenerate (CLOSE) $\times2$',
        'color': paultol.muted.colors[3],
        'linestyle': '-',
        'offset': 2,
    }
    # Comparison DTDs
    plateau = {
        'func': dtds.plateau(width=0.3, slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=0.3$ Gyr) $\times2$',
        'color': paultol.muted.colors[7],
        'linestyle': '--',
        'offset': 2,
    }
    plateau_long = {
        'func': dtds.plateau(width=1., slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=1$ Gyr)',
        'color': paultol.muted.colors[6],
        'linestyle': '--',
        'offset': 1,
    }
    exp = {
        'func': dtds.exponential(timescale=1.5, tmin=MIN_RIA_DELAY),
        'label': r'Exponential ($\tau=1.5$ Gyr) $\times6$',
        'color': paultol.muted.colors[0],
        'linestyle': '--',
        'offset': 4,
    }
    distlist = [sd, exp, dd_close, plateau, dd_wide, plateau_long]


if __name__ == '__main__':
    main()
