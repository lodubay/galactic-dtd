"""
Plot a selection of models for the Type Ia supernova delay time distributions 
(DTD) as a function of time.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import paths
from multizone.src import dtds
from colormaps import paultol
from _globals import MIN_RIA_DELAY, ONE_COLUMN_WIDTH

def main(style='paper'):
    plt.style.use(paths.styles / f'{style}.mplstyle')
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)
    # get default line width
    default_line_width = plt.rcParams['lines.linewidth']
    fig, ax = setup_axes()
    times = np.arange(0.04, 13.2, 0.001)
    distributions = [styles.prompt, styles.plaw, styles.exp, 
                     styles.plateau_long, styles.triple]
    for dtd in distributions:
        func = dtd['func']
        yvals = np.array([func(t) for t in times])
        ax.plot(times, yvals / func(1.), 
                label=dtd['label'], c=dtd['color'], ls=dtd['line'], 
                lw=default_line_width * dtd['lwmod'])
        # indicate median delay times
        cdf = np.cumsum(yvals / np.sum(yvals))
        med_idx = np.where(cdf >= 0.5)[0][0]
        med = times[med_idx]
        ax.scatter(med, 5e-3, c=dtd['color'], s=10, marker='o')
    ax.set_ylim((3e-3, 3e2))
    # Plot the SDSS-II DTD recovered by Maoz et al. (2012), MNRAS 426, 3282
    # (see their Table 2)
    tbins = np.array([0.04, 0.42, 2.4, 14])
    times = (tbins[:-1] + tbins[1:])/2
    tbin_widths = (tbins[1:] - tbins[:-1])/2
    sn_rates = np.array([140, 25.1, 1.83])
    sn_rate_errors = np.array([30, 6.3, 0.42])
    # Linear regression
    A = np.vstack([np.log10(times), np.ones(len(times))]).T
    m, c = np.linalg.lstsq(A, np.log10(sn_rates), rcond=None)[0]
    # Scale so best-fit line passes through (1, 1)
    norm = 1 / 10**c
    ax.errorbar(times, sn_rates*norm, xerr=tbin_widths, yerr=sn_rate_errors*norm,
                color='k', linestyle='none', capsize=1, elinewidth=0.5,
                capthick=0.5, marker='s', markersize=2, 
                label='Maoz et al. (2012)')
    
    ax.legend(frameon=False, loc='upper right', handlelength=2)
    fig.savefig(paths.figures / 'delay_time_distributions')
    plt.close()


class styles:
    """Plot styling for different DTD models."""
    plaw = {
        'func': dtds.powerlaw(slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Power-law ($\alpha=-1.1$)',
        'color': paultol.bright.colors[5], # purple
        'line': '-',
        'lwmod': 1.2,
    }
    plateau = {
        'func': dtds.plateau(width=0.3, slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=0.3$ Gyr)',
        'color': paultol.bright.colors[2], # green
        'line': '--',
        'lwmod': 1.2,
    }
    exp = {
        'func': dtds.exponential(timescale=1.5, tmin=MIN_RIA_DELAY),
        'label': r'Exponential ($\tau=1.5$ Gyr)',
        'color': paultol.bright.colors[0], # blue
        'line': (5, (10, 3)), # long dashed with offset
        'lwmod': 1.2,
    }
    prompt = {
        'func': dtds.prompt(peak=0.05, stdev=0.015, timescale=3, tmin=MIN_RIA_DELAY),
        'label': 'Two-population',
        'color': paultol.bright.colors[1], # red
        'line': ':',
        'lwmod': 1.8,
    }
    triple = {
        'func': dtds.triple(tmin=MIN_RIA_DELAY),
        'label': r'Triple-system',
        'color': paultol.bright.colors[3], # yellow
        'line': '-.',
        'lwmod': 1.2,
    }
    # Additional DTDs
    exp_long = {
        'func': dtds.exponential(timescale=3, tmin=MIN_RIA_DELAY),
        'label': r'Exponential ($\tau=3$ Gyr)',
        'color': paultol.bright.colors[0], # blue
        'line': '--',
        'lwmod': 1.2,
    }
    plateau_long = {
        'func': dtds.plateau(width=1., slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': r'Plateau ($W=1$ Gyr)',
        'color': paultol.bright.colors[2], # green
        'line': '--',
        'lwmod': 1.2,
    }
    plaw_steep = {
        'func': dtds.powerlaw(slope=-1.4, tmin=MIN_RIA_DELAY),
        'label': r'Power-law ($\alpha=-1.4$)',
        'color': paultol.bright.colors[5], # purple
        'line': '--',
        'lwmod': 1.2,
    }


def setup_axes(width=ONE_COLUMN_WIDTH):
    fig, ax = plt.subplots(figsize=(width, width))
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.1, top=0.98)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.set_xlabel('Time after star formation [Gyr]')
    ax.set_ylabel('Relative SN Ia rate')
    return fig, ax
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='delay_time_distributions.py',
        description='Plot a selection of the delay time distribution models ' +
                    'as a function of time',
        )
    parser.add_argument('-s', '--style', 
                        choices=['paper', 'poster'],
                        default='paper', 
                        help='Plot style to use (default: paper)')
    args = parser.parse_args()
    main(**vars(args))
