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

SD_FACTOR = 1.2

def main():
    fig, ax = plt.subplots()
    tarr = np.arange(0.04, END_TIME, 0.001)
    # WIDE DD
    dd_wide = dtds.greggio05_double(scheme='wide', dt=1e-3, nsamples=50)
    # # tarr = dd_wide.times
    dd_wide_arr = np.array([dd_wide(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_wide_arr, label='DD WIDE')
    # # Approximation to WIDE DD
    dd_wide_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_wide_arr)
    ax.plot(tarr * 1e9, [dd_wide_fit(t) for t in tarr], ls='--',
            label='Approximate function for wide DD')
    # # Fit to WIDE DD
    slope1 = linregress(np.log10(tarr[(tarr >= 0.2) & (tarr < 1)]),
                        np.log10(dd_wide_arr[(tarr >= 0.2) & (tarr < 1)]))[0]
    print(slope1)
    slope2 = linregress(np.log10(tarr[tarr >= 2]),
                        np.log10(dd_wide_arr[tarr >= 2]))[0]
    print(slope2)
    bp = dtds.powerlaw_broken(slope1=slope1, slope2=slope2, tsplit=1)
    ax.plot(tarr * 1e9, [bp(t) for t in tarr],
            label='Broken powerlaw fit to close DD')

    # CLOSE DD
    dd_close = dtds.greggio05_double(scheme='close', dt=1e-3, nsamples=100)
    dd_close_arr = np.array([dd_close(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_close_arr, label='DD CLOSE')
    # Approximation to CLOSE DD
    dd_close_fit = dtds.greggio05_approximate.fit_to_data(tarr, dd_close_arr)
    print('Slope 1 = %s\nSlope 2 = %s' % (dd_close_fit.slope1, dd_close_fit.slope2))
    print('Rise strength = %s, timescale = %s' % (dd_close_fit.rise.strength, dd_close_fit.rise.timescale))
    print('Tail strength = %s, timescale = %s' % (dd_close_fit.tail.strength, dd_close_fit.tail.timescale))
    # dd_close_fit = dtds.greggio05_approximate.from_defaults('close')
    ax.plot(tarr * 1e9, [dd_close_fit(t) for t in tarr], ls='--',
            label='Approximate function for wide DD')
    # popt, pcov = curve_fit(dd_approx_func, tarr, dd_close_arr,
    #                         p0=(1, 0.1, 0, -1, 1, 0.1, 1e-9))
    # print(popt)
    # ax.plot(tarr * 1e9, dd_approx_func(tarr, *popt), ls='--',
    #         label='Approximate function for close DD')
    # Fit to CLOSE DD
    # popt, pcov = curve_fit(broken_powerlaw, tarr[tarr > 0.2],
    #                        dd_close_arr[tarr > 0.2],
    #                        p0=(1e-9, 0, 1, -1))
    # print(popt)
    # ax.plot(tarr * 1e9, broken_powerlaw(tarr, *popt),
    #         label='Broken powerlaw fit to close DD')
    # popt, pcov = curve_fit(powerlaw, tarr[(tarr >= 0.2) & (tarr < 1)],
    #                        dd_close_arr[(tarr >= 0.2) & (tarr < 1)])
    # slope1 = popt[0]
    # print(slope1)
    # popt, pcov = curve_fit(powerlaw, tarr[tarr >= 2],
    #                        dd_close_arr[tarr >= 2])
    # slope2 = popt[0]
    # print(slope2)
    slope1 = linregress(np.log10(tarr[(tarr >= 0.2) & (tarr < 1)]),
                        np.log10(dd_close_arr[(tarr >= 0.2) & (tarr < 1)]))[0]
    print(slope1)
    slope2 = linregress(np.log10(tarr[tarr >= 2]),
                        np.log10(dd_close_arr[tarr >= 2]))[0]
    print(slope2)
    bp = dtds.powerlaw_broken(slope1=slope1, slope2=slope2, tsplit=1)
    ax.plot(tarr * 1e9, [bp(t) for t in tarr],
            label='Broken powerlaw fit to close DD')

    # SD
    sd = dtds.greggio05_single()
    sd_arr = np.array([sd(t) * SD_FACTOR for t in tarr])
    ax.plot(tarr * 1e9, sd_arr, label=r'SD $\times$ factor')

    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.04e9, 16e9))
    ax.set_ylabel(r'$f_{\rm{Ia}}$')
    ax.set_yscale('log')
    ax.set_ylim((3e-12, 1e-9))
    ax.legend()
    fig.savefig(paths.figures / 'greggio05.png', dpi=300)

if __name__ == '__main__':
    main()
