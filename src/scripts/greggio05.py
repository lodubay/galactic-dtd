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

SD_FACTOR = 2
WIDE_COLORS = [paultol.muted.colors[1],
               paultol.muted.colors[4],
               paultol.muted.colors[6]]
CLOSE_COLORS = [paultol.muted.colors[3],
                paultol.muted.colors[2],
                paultol.muted.colors[7]]

def main(scheme, dt=1e-3, nsamples=100):
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.25), tight_layout=True,
                            sharex=True)
    tarr = np.arange(0.04, END_TIME, dt)

    if scheme == 'wide':
        colors = WIDE_COLORS
    else:
        colors = CLOSE_COLORS

    # Left panel: RIa
    ax = axs[0]

    # Single-degenerate DTD for comparison
    print('Plotting single-degenerate DTD...')
    sd = dtds.greggio05_single()
    sd_arr = np.array([sd(t) for t in tarr])
    ax.plot(tarr * 1e9, SD_FACTOR * sd_arr,
            color=paultol.muted.colors[0],
            label=r'SD $\times$ %s' % SD_FACTOR)
    print('Done!')

    # The original Greggio 2005 DD DTD (takes a while if nsamples is large!)
    dd_func = dtds.greggio05_double(scheme, dt=dt, nsamples=nsamples)
    dd_arr = np.array([dd_func(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_arr,
            color=colors[0], linestyle='-',
            label='DD %s' % scheme.upper())

    # Best fit approximate function
    dd_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_arr)
    fit_arr = np.array([dd_fit(t) for t in tarr])
    ax.plot(tarr * 1e9, fit_arr,
            color=colors[1], linestyle='--',
            label='Best-fit approximation')

    # Broken power-law linear regression
    print('Computing broken power-law linear regression...')
    slope1 = linregress(np.log10(tarr[(tarr >= 0.2) & (tarr < 1)]),
                        np.log10(dd_arr[(tarr >= 0.2) & (tarr < 1)]))[0]
    slope2 = linregress(np.log10(tarr[tarr >= 2]),
                        np.log10(dd_arr[tarr >= 2]))[0]
    bp = dtds.powerlaw_broken(slope1=slope1, slope2=slope2, tsplit=1)
    bp_arr = np.array([bp(t) for t in tarr])
    ax.plot(tarr * 1e9, bp_arr,
            color=colors[2], linestyle=':',
            label='Broken power-law fit')
    # Add text with power-law slopes
    t1 = 0.7
    ax.text(t1 * 1e9, 0.6 * bp(t1), rf'$\alpha_1={slope1:.02f}$', ha='right',
            fontsize=8)
    t2 = 3
    ax.text(t2 * 1e9, 0.5 * bp(t2), rf'$\alpha_2={slope2:.02f}$', ha='right',
            fontsize=8)
    print('Done!')

    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.03e9, 16e9))
    ax.set_ylabel(r'$R_{\rm{Ia}}$')
    ax.set_yscale('log')
    ax.set_ylim((3e-12, 3e-9))
    # ax.legend(fontsize=8, frameon=False, handlelength=1.5)
    handles, labels = ax.get_legend_handles_labels()

    # Right panel: CDF
    ax = axs[1]
    ax.plot(tarr * 1e9, 1e9 * dt * np.cumsum(sd_arr), c=paultol.muted.colors[0],
            ls='-')
    ax.plot(tarr * 1e9, 1e9 * dt * np.cumsum(dd_arr), c=colors[0], ls='-')
    ax.plot(tarr * 1e9, 1e9 * dt * np.cumsum(fit_arr), c=colors[1], ls='--')
    ax.plot(tarr * 1e9, 1e9 * dt * np.cumsum(bp_arr), c=colors[2], ls=':')
    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('CDF')
    ax.set_yscale('log')
    ax.set_ylim((3e-4, 2))
    ax.legend(handles, labels, loc='lower right', fontsize=8, frameon=False,
              handlelength=1.5)

    fig.savefig(paths.figures / ('greggio05_%s.png' % scheme), dpi=300)

if __name__ == '__main__':
    main('close')
