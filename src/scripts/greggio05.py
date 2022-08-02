# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:06:58 2022

@author: dubay.11
"""

import sys
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from migration.src._globals import END_TIME
from colormaps import paultol

SD_FACTOR = 1
WIDE_COLORS = [paultol.muted.colors[1],
               paultol.muted.colors[4],
               paultol.muted.colors[6]]
CLOSE_COLORS = [paultol.muted.colors[3],
                paultol.muted.colors[2],
                paultol.muted.colors[7]]

def main(scheme, dt=1e-4, nsamples=200, verbose=False, format='pdf'):
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    tarr = np.arange(0.04, END_TIME, dt)

    if scheme == 'wide':
        colors = WIDE_COLORS
    else:
        colors = CLOSE_COLORS

    # Single-degenerate DTD for comparison
    if verbose:
        print('Plotting single-degenerate DTD...')
    sd = dtds.greggio05_single()
    sd_arr = np.array([sd(t) for t in tarr])
    ax.plot(tarr * 1e9, SD_FACTOR * sd_arr,
            color=paultol.muted.colors[0], linestyle='-', linewidth=1,
            label=r'SD $\times$ %s' % SD_FACTOR)
    if verbose:
        print('Done!')

    # The original Greggio 2005 DD DTD (takes a while if nsamples is large!)
    dd_func = dtds.greggio05_double(scheme, dt=dt, nsamples=nsamples,
                                    progress=verbose)
    dd_arr = np.array([dd_func(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_arr,
            color=colors[0], linestyle='-', linewidth=1,
            label='DD %s' % scheme.upper())

    # Best fit approximate function
    if verbose:
        print('Fitting approximate model to Greggio 2005...')
    dd_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_arr)
    fit_arr = np.array([dd_fit(t) for t in tarr])
    ax.plot(tarr * 1e9, fit_arr,
            color=colors[1], linestyle='--', linewidth=1,
            label='Best-fit approximation')
    if verbose:
        print('Done!')

    # Broken power-law linear regression
    if verbose:
        print('Computing broken power-law linear regression...')
    slope1 = linregress(np.log10(tarr[(tarr >= 0.2) & (tarr < 1)]),
                        np.log10(dd_arr[(tarr >= 0.2) & (tarr < 1)]))[0]
    slope2 = linregress(np.log10(tarr[tarr >= 2]),
                        np.log10(dd_arr[tarr >= 2]))[0]
    bp = dtds.powerlaw_broken(slope1=slope1, slope2=slope2, tsplit=1)
    bp_arr = np.array([bp(t) for t in tarr])
    ax.plot(tarr * 1e9, bp_arr,
            color=colors[2], linestyle=':', linewidth=1,
            label='Broken power-law fit')
    # Add text with power-law slopes
    t1 = 0.5
    ax.text(t1 * 1e9, 0.7 * bp(t1), rf'$\alpha_1={slope1:.02f}$', ha='right',
            fontsize=7)
    t2 = 3
    ax.text(t2 * 1e9, 1.3 * bp(t2), rf'$\alpha_2={slope2:.02f}$', ha='left',
            fontsize=7)
    if verbose:
        print('Done!')

    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.03e9, 16e9))
    ax.set_ylabel(r'Normalized SN Ia Rate [M$^{-1}_{\odot}$ yr$^{-1}$]')
    ax.set_yscale('log')
    ax.set_ylim((3e-12, 3e-9))
    ax.legend(loc='upper right', fontsize=7, frameon=False, handlelength=1.5)

    fig.savefig(paths.figures / ('greggio05_%s.%s' % (scheme, format)), dpi=300)

if __name__ == '__main__':
    main('close', format='png', verbose=True, dt=1e-3, nsamples=100)
