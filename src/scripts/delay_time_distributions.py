"""
Plot the Type Ia supernova delay time distributions (DTDs) as a function of time.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import paths
from multizone.src import dtds
from colormaps import paultol
from _globals import MIN_RIA_DELAY
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig, ax= setup_axes()
    times = [t*0.001 for t in range(40, 13200)]
    distributions = [styles.prompt, styles.plaw, styles.plateau, 
                      styles.exp, styles.triple]
    for dtd in distributions:
        func = dtd['func']
        ax.plot(times, [func(t) / func(1) for t in times], 
                label=dtd['label'], c=dtd['color'], ls=dtd['line'], lw=1)
    ax.set_ylim((3e-3, 3e2))
    ax.legend(frameon=False, loc='upper right', fontsize=8, handlelength=1.25)
    fig.savefig(paths.figures / 'delay_time_distributions.pdf')
    plt.close()


class styles:
    """Plot styling for different DTD models."""
    plaw = {
        'func': dtds.powerlaw(slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Power-law ($\alpha=-1.1$)',
        'color': paultol.bright.colors[0],
        'line': '-',
    }
    plateau = {
        'func': dtds.plateau(width=0.3, slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=0.3$ Gyr)',
        'color': paultol.bright.colors[4],
        'line': '--',
    }
    exp = {
        'func': dtds.exponential(timescale=1.5, tmin=MIN_RIA_DELAY),
        'label': r'Exponential ($\tau=1.5$ Gyr)',
        'color': paultol.bright.colors[5],
        'line': '-',
    }
    prompt = {
        'func': dtds.prompt(peak=0.05, stdev=0.015, timescale=3, tmin=MIN_RIA_DELAY),
        'label': 'Two-population',
        'color': paultol.bright.colors[1],
        'line': ':',
    }
    triple = {
        'func': dtds.triple(tmin=MIN_RIA_DELAY),
        'label': r'Triple-system',
        'color': paultol.bright.colors[2],
        'line': '-.'
    }
    # Additional DTDs
    exp_long = {
        'func': dtds.exponential(timescale=3, tmin=MIN_RIA_DELAY),
        'label': r'Exponential ($\tau=3$ Gyr)',
        'color': paultol.bright.colors[5],
        'line': '--',
    }
    plateau_long = {
        'func': dtds.plateau(width=1., slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=1$ Gyr)',
        'color': paultol.bright.colors[4],
        'line': ':',
    }
    plaw_steep = {
        'func': dtds.powerlaw(slope=-1.4, tmin=MIN_RIA_DELAY),
        'label': r'Power-law ($\alpha=-1.4$)',
        'color': paultol.bright.colors[0],
        'line': '--',
    }


def setup_axes():
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.set_xlabel('Time after star formation [Gyr]')
    ax.set_ylabel('Relative SN Ia rate')
    return fig, ax
    

if __name__ == '__main__':
    main()
