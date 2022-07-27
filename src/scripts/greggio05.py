# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:06:58 2022

@author: dubay.11
"""

import sys
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from migration.src._globals import END_TIME
from colormaps import paultol

def main():
    fig, ax = plt.subplots()
    # tarr = np.logspace(np.log10(0.04), np.log10(END_TIME), num=200, base=10)
    # plaw_slope = -0.8
    # broken_powerlaw = dtds.powerlaw_broken(tsplit=1, slope2=plaw_slope)
    # ax.plot(tarr * 1e9, [broken_powerlaw(t, slope=plaw_slope) for t in tarr],
    #         color='k',
    #         label=r'Broken powerlaw ($\alpha=%s$)' % plaw_slope)
    # for scheme, color in zip(['wide', 'close'], ['b', 'g']):
    #     for slope, ls in zip([-1.44, -0.7], ['-', '--']):
    #         print(' '.join((scheme, str(slope))))
    #         dtd = dtds.greggio05(scheme=scheme, sd_slope=slope)
    #         RIa = [dtd(t) for t in tarr]
    #         ax.plot(tarr * 1e9, RIa, color=color, ls=ls,
    #                  label=r'%s scheme, $|\dot m_2|\propto\tau^{%s}$' % (
    #                      scheme.upper(), slope))

    dd_wide = dtds.greggio05_double(scheme='wide', dt=1e-3, nsamples=200,
                                    sd_slope=-0.44, efficiency=1)
    tarr = dd_wide.times
    ax.plot(tarr * 1e9, [2*dd_wide(t) for t in tarr], label='DD WIDE')
    dd_close = dtds.greggio05_double(scheme='close', dt=1e-3, nsamples=200,
                                     sd_slope=-0.44, efficiency=1)
    ax.plot(tarr * 1e9, [2*dd_close(t) for t in tarr], label='DD CLOSE')
    sd = dtds.greggio05_single(m2_slope=-0.44, efficiency=1)
    ax.plot(tarr * 1e9, [27*sd(t) for t in tarr], label='SD')
    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.04e9, 16e9))
    ax.set_ylabel(r'$f_{\rm{Ia}}$')
    ax.set_yscale('log')
    ax.set_ylim((1e-3, 10))
    fig.legend()
    fig.savefig(paths.figures / 'greggio05.png', dpi=300)

def broken_powerlaw(time, slope=-0.8):
    if time >= 0.05 and time < 1:
        return 0.5
    elif time >= 1:
        return 0.2 * time ** slope
    else:
        return 0

if __name__ == '__main__':
    main()
