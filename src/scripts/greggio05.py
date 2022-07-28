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
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from migration.src._globals import END_TIME
from colormaps import paultol

def main():
    fig, ax = plt.subplots()
    tarr = np.arange(0.04, END_TIME, 0.001)
    # WIDE DD
    # dd_wide = dtds.greggio05_double(scheme='wide', dt=1e-3, nsamples=100)
    # # tarr = dd_wide.times
    # dd_wide_arr = np.array([dd_wide(t) for t in tarr])
    # ax.plot(tarr * 1e9, dd_wide_arr, label='DD WIDE')
    # # Approximation to WIDE DD
    # popt, pcov = curve_fit(dd_approx_func, tarr, dd_wide_arr,
    #                        p0=(1, 0.1, 0, -1, 1, 0.1, 1e-9))
    # print(popt)
    # ax.plot(tarr * 1e9, dd_approx_func(tarr, *popt),
    #         label='Approximate function for wide DD')
    # # Fit to WIDE DD
    # popt, pcov = curve_fit(broken_powerlaw, tarr[tarr > 0.2],
    #                        dd_wide_arr[tarr > 0.2],
    #                        p0=(1e-9, 0, 1, -1))
    # print(popt)
    # ax.plot(tarr * 1e9, broken_powerlaw(tarr, *popt),
    #         label='Broken powerlaw fit to wide DD')

    # CLOSE DD
    dd_close = dtds.greggio05_double(scheme='close', dt=1e-3, nsamples=100)
    dd_close_arr = np.array([dd_close(t) for t in tarr])
    ax.plot(tarr * 1e9, dd_close_arr, label='DD CLOSE')
    # Approximation to CLOSE DD
    popt, pcov = curve_fit(dd_approx_func, tarr, dd_close_arr,
                            p0=(1, 0.1, 0, -1, 1, 0.1, 1e-9))
    print(popt)
    ax.plot(tarr * 1e9, dd_approx_func(tarr, *popt),
            label='Approximate function for close DD')
    # Fit to CLOSE DD
    popt, pcov = curve_fit(broken_powerlaw, tarr, dd_close_arr,
                            p0=(1e-9, 0, 1, -1))
    print(popt)
    ax.plot(tarr * 1e9, broken_powerlaw(tarr, *popt),
            label='Broken powerlaw fit to close DD')

    # SD
    # sd = dtds.greggio05_single()
    # sd_arr = np.array([sd(t) for t in tarr])
    # ax.plot(tarr * 1e9, sd_arr, label='SD')

    ax.set_xlabel('Time [yr]')
    ax.set_xscale('log')
    ax.set_xlim((0.04e9, 16e9))
    ax.set_ylabel(r'$f_{\rm{Ia}}$')
    ax.set_yscale('log')
    # ax.set_ylim((1e-3, 10))
    ax.legend()
    fig.savefig(paths.figures / 'greggio05.png', dpi=300)

def broken_powerlaw(time, norm, slope1, tsplit, slope2):
    part1 = time[time <= tsplit] ** slope1
    part2 = tsplit ** (slope1 - slope2) * time[time > tsplit] ** slope2
    return norm * np.concatenate((part1, part2))

def dd_approx_func(time, rise_strength, rise_timescale, slope1, slope2,
                tail_strength, tail_timescale, norm):
    approx_class = dd_approx(rise_strength=rise_strength,
                             rise_timescale=rise_timescale,
                             slope1=slope1, slope2=slope2,
                             tail_strength=tail_strength,
                             tail_timescale=tail_timescale,
                             norm=norm)
    if isinstance(time, np.ndarray):
        return np.array([approx_class(t) for t in time])
    else:
        return approx_class(time)

class dd_approx:
    def __init__(self, rise_strength=3, rise_timescale=0.05, slope1=-0.3,
                 tsplit=1, slope2=-1.3, tail_strength=1, tail_timescale=0.1,
                 norm=3e-10):
        self.rise_strength = rise_strength
        self.rise_timescale = rise_timescale
        self.slope1= slope1
        self.tsplit = tsplit
        self.slope2 = slope2
        self.tail_strength = tail_strength
        self.tail_timescale = tail_timescale
        self.norm = norm

    def __call__(self, time):
        if time <= self.tsplit:
            return self.norm * self.part1(time)
        else:
            return self.norm * \
                (self.part1(self.tsplit) / self.part2(self.tsplit)) * \
                self.part2(time)

    def part1(self, time):
        return self.exp_rise(time) * time ** self.slope1

    def exp_rise(self, time):
        return (1 - self.rise_strength * np.exp(-time/self.rise_timescale))

    def part2(self, time):
        return self.exp_tail(time) * time ** self.slope2

    def exp_tail(self, time):
        return (1 + self.tail_strength * np.exp(-(time-self.tsplit)/self.tail_timescale))

def sd_approx_func(time, rise_strength, rise_timescale, norm):
    approx_class = sd_approx(rise_strength=rise_strength,
                             rise_timescale=rise_timescale,
                             timescale=1.5, norm=norm)
    if isinstance(time, np.ndarray):
        return np.array([approx_class(t) for t in time])
    else:
        return approx_class(time)

class sd_approx:
    def __init__(self, rise_strength=3, rise_timescale=0.05, timescale=1.5,
                 norm=1e-9):
        self.rise_strength = rise_strength
        self.rise_timescale = rise_timescale
        self.timescale = timescale
        self.norm = norm

    def __call__(self, time):
        return (1 - self.rise_strength * np.exp(-time/self.rise_timescale)) * \
            np.exp(-time/self.timescale) * self.norm

# def broken_powerlaw(time, tsplit, slope1, slope2):
#     bp_class = dtds.powerlaw_broken(tsplit=tsplit, slope1=slope1, slope2=slope2)
#     if isinstance(time, np.ndarray):
#         return np.array([bp_class(t) for t in time])
#     elif isinstance(time, float):
#         return bp_class(time)
#     else:
#         raise TypeError('Time must be a float or array of floats')

if __name__ == '__main__':
    main()
