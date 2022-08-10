# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:06:58 2022

@author: dubay.11
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from migration.src._globals import END_TIME
from colormaps import paultol

DT = 1e-4
NSAMPLES = 200

def main(dt=DT, nsamples=NSAMPLES, verbose=True):
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    tarr = np.arange(0.04, END_TIME, dt)

    # Single-degenerate DTD
    if verbose:
        print('Plotting single-degenerate DTD...')
    sd = dtds.greggio05_single()
    sd_arr = np.array([sd(t) for t in tarr])
    ax.plot(tarr * 1e9, sd_arr,
            color=paultol.muted.colors[0], linestyle='-', linewidth=1,
            label='Single Degenerate' )
    if verbose:
        print('Done!')
        
    scheme = 'close'

    # The original Greggio 2005 DD DTD (takes a while if nsamples is large!)
    dd_func = dtds.greggio05_double(scheme, dt=dt, nsamples=nsamples,
                                    progress=verbose)
    dd_arr = np.array([dd_func(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_arr,
            color=paultol.muted.colors[3], linestyle='-', linewidth=1,
            label='Double Degenerate %s' % scheme.upper())

    # Best fit approximate function
    if verbose:
        print('Fitting approximate model to DD %s...' % scheme.upper())
    dd_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_arr)
    fit_arr = np.array([dd_fit(t) for t in tarr])
    ax.plot(tarr * 1e9, fit_arr,
            color=paultol.muted.colors[2], linestyle='--', linewidth=1,
            label='DD %s Best-Fit Approximation' % scheme.upper())
    if verbose:
        print('Done!')
        
    scheme = 'wide'

    # The original Greggio 2005 DD DTD (takes a while if nsamples is large!)
    dd_func = dtds.greggio05_double(scheme, dt=dt, nsamples=nsamples,
                                    progress=verbose)
    dd_arr = np.array([dd_func(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_arr,
            color=paultol.muted.colors[1], linestyle='-', linewidth=1,
            label='Double Degenerate %s' % scheme.upper())

    # Best fit approximate function
    if verbose:
        print('Fitting approximate model to DD %s...' % scheme.upper())
    dd_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_arr)
    fit_arr = np.array([dd_fit(t) for t in tarr])
    ax.plot(tarr * 1e9, fit_arr,
            color=paultol.muted.colors[4], linestyle='--', linewidth=1,
            label='DD %s Best-Fit Approximation' % scheme.upper())
    if verbose:
        print('Done!')

    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.03e9, 16e9))
    ax.set_ylabel(r'Normalized SN Ia Rate [M$^{-1}_{\odot}$ yr$^{-1}$]')
    ax.set_yscale('log')
    ax.set_ylim((2e-12, 2e-9))
    ax.legend(loc='lower center', fontsize=7, handlelength=1.2, frameon=False)

    fig.savefig(paths.figures / 'greggio05.pdf', dpi=300)

if __name__ == '__main__':
    main()
