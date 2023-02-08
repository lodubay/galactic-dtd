"""
Plot analytical DTDs from Greggio (2005) alongside closest simpler models.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from migration.src._globals import END_TIME
from delay_time_distributions import setup_axes
from colormaps import paultol

DT = 1e-4 # plotting timestep in Gyr
NSAMPLES = 200 # sets the integration precision; higher means longer runtime
DELAY = 0.04 # minimum SN Ia delay time in Gyr

def main(dt=DT, nsamples=NSAMPLES, verbose=True):
    
    class styles:
        """Plot styling for different DTD models."""
        # Greggio (2005) DTDs
        sd = {
            'func': dtds.greggio05_single(),
            'label': 'Single Degenerate',
            'color': paultol.muted.colors[5],
            'line': '-'
        }
        dd_wide = {
            # 'func': dtds.greggio05_double('wide', dt=dt, nsamples=nsamples,
            #                               progress=verbose),
            'func': dtds.greggio05_approximate.from_defaults('wide'),
            'label': 'Double Degenerate (WIDE)',
            'color': paultol.muted.colors[1],
            'line': '-',
        }
        dd_close = {
            # 'func': dtds.greggio05_double('close', dt=dt, nsamples=nsamples,
            #                               progress=verbose),
            'func': dtds.greggio05_approximate.from_defaults('close'),
            'label': 'Double Degenerate (CLOSE)',
            'color': paultol.muted.colors[3],
            'line': '-',
        }
        # Comparison DTDs
        plateau = {
            'func': dtds.plateau(width=0.3, slope=-1.1, tmin=DELAY),
            'label': r'Plateau ($W=0.3$ Gyr)',
            'color': paultol.muted.colors[7],
            'line': '--',
        }
        plateau_long = {
            'func': dtds.plateau(width=1., slope=-1.1, tmin=DELAY),
            'label': r'Plateau ($W=1$ Gyr)',
            'color': paultol.muted.colors[6],
            'line': '--',
        }
        exp = {
            'func': dtds.exponential(timescale=1.5, tmin=DELAY),
            'label': r'Exponential ($\tau=1.5$ Gyr)',
            'color': paultol.muted.colors[0],
            'line': '--',
        }
    
    fig, ax = setup_axes()
    tarr = np.arange(0.04, END_TIME, dt)
    distributions = [styles.sd, styles.dd_close, styles.dd_wide, styles.plateau,
                     styles.exp, styles.plateau_long]
    
    for dtd in distributions:
        if verbose:
            print('Plotting', dtd['label'], 'DTD')
        func = dtd['func']
        ax.plot(tarr, np.array([func(t) for t in tarr]), 
                label=dtd['label'], c=dtd['color'], ls=dtd['line'], lw=1)

    ax.set_ylim((3e-2, 3e1))
    ax.legend(loc='upper right', fontsize=8, handlelength=1.2, frameon=False)

    fig.savefig(paths.figures / 'dtd_greggio05.png', dpi=300)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='dtd_greggio05.py',
        description='Plot analytical DTDs from Greggio (2005).'
    )
    parser.add_argument('--dt', metavar='DT', type=float, default=DT,
                        help='Plotting timestep in Gyr')
    parser.add_argument('--nsamples', '-n', metavar='NSAMPLES', type=int,
                        default=NSAMPLES,
                        help='Numerical integration precision (higher ' + \
                             'takes longer to plot)')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    main(**vars(args))
