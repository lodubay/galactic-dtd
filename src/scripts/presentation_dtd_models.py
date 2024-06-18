"""
Simplified plot comparing the different DTD models for presentation purposes.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import paths
from multizone.src import dtds
from colormaps import paultol
from delay_time_distributions import setup_axes
from _globals import MIN_RIA_DELAY

def main():
    plt.style.use(paths.styles / 'presentation.mplstyle')
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)
    # get default line width
    default_line_width = plt.rcParams['lines.linewidth']
    fig, ax = setup_axes(width=6)
    times = np.arange(0.04, 13.2, 0.001)
    distributions = [styles.plaw, styles.exp, styles.plateau_long]
    for dtd in distributions:
        func = dtd['func']
        yvals = np.array([func(t) for t in times])
        ax.plot(times, yvals / func(1.), 
                label=dtd['label'], c=dtd['color'], ls=dtd['line'], 
                lw=default_line_width)
        # indicate median delay times
        cdf = np.cumsum(yvals / np.sum(yvals))
        med_idx = np.where(cdf >= 0.5)[0][0]
        med = times[med_idx]
        ax.scatter(med, 0.013, c=dtd['color'], s=40, marker='o')
    # Label median delay times
    ax.text(1, 0.017, 'Median delay times', ha='center')
    ax.set_ylim((1e-2, 50))
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
                color='k', linestyle='none', capsize=2, 
                elinewidth=0.5 * default_line_width,
                capthick=0.5 * default_line_width, marker='s', markersize=5, 
                label='SN survey recovery\n(Maoz et al. 2012)')
    
    # Replace y-axis label
    ax.set_ylabel('Relative supernova rate', labelpad=-6)
    ax.legend(frameon=False, loc='upper right', handlelength=1.2)
    fig.savefig(paths.extra / 'presentation' / 'dtd_models')
    plt.close()


class styles:
    """Plot styling for different DTD models."""
    plaw = {
        'func': dtds.powerlaw(slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': 'Observational',
        'color': paultol.bright.colors[5], # purple
        'line': '-',
    }
    exp = {
        'func': dtds.exponential(timescale=1.5, tmin=MIN_RIA_DELAY),
        'label': 'Theoretical WD+star',
        'color': paultol.bright.colors[0], # blue
        'line': '--',
    }
    plateau_long = {
        'func': dtds.plateau(width=1., slope=-1.1, tmin=MIN_RIA_DELAY),
        'label': 'Theoretical WD+WD',
        'color': paultol.bright.colors[2], # green
        'line': '-.',
    }
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='presentation_dtd_models.py',
        description='Presentation-ready version of a plot comparing DTD models',
        )
    args = parser.parse_args()
    main(**vars(args))