"""
This script plots one-zone tracks for the four star formation history (SFH)
models in the paper.
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from colormaps import paultol
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)
from multizone.src.yields import J21
from multizone.src import models, dtds
from track_and_mdf import setup_axes, plot_vice_onezone
from _globals import END_TIME, ONEZONE_DEFAULTS

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    
    output_dir = paths.data / 'onezone' / 'sfh'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    rgal = 8.
    dr = 0.1
    ONEZONE_DEFAULTS['RIa'] = dtds.exponential(timescale=1.5)
    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    
    # Inside-out SFH
    sz = vice.singlezone(name=str(output_dir / 'insideout'),
                         func=models.insideout(rgal, dr=dr, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    
    # Late-burst SFH
    sz = vice.singlezone(name=str(output_dir / 'lateburst'),
                         func=models.lateburst(rgal, dr=dr, dt=dt), 
                         mode='sfr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    
    # Two-infall SFH
    sz = vice.singlezone(name=str(output_dir / 'twoinfall'),
                         func=models.twoinfall(rgal, dr=dr, dt=dt), 
                         mode='ifr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    
    # Early-burst SFH (no area / Mgas dependence)
    tau_star = models.earlyburst_tau_star(1)
    ONEZONE_DEFAULTS['tau_star'] = tau_star.time_dependence
    sz = vice.singlezone(name=str(output_dir / 'earlyburst'),
                         func=models.earlyburst_ifr(rgal, dr=dr, dt=dt), 
                         mode='ifr',
                         **ONEZONE_DEFAULTS)
    sz.run(simtime, overwrite=True)
    
    # Plot
    fig, axs = setup_axes()
    plot_vice_onezone(str(output_dir / 'insideout'), 
                      fig=fig, axs=axs, label='Inside-out',
                      style_kw={'zorder': 1, 'linestyle': '-'})
    plot_vice_onezone(str(output_dir / 'lateburst'), 
                      fig=fig, axs=axs, label='Late-burst',
                      style_kw={'zorder': 2, 'linestyle': '--'})
    plot_vice_onezone(str(output_dir / 'earlyburst'), 
                      fig=fig, axs=axs, label='Early-burst',
                      style_kw={'zorder': 3, 'linestyle': ':'})
    plot_vice_onezone(str(output_dir / 'twoinfall'), 
                      fig=fig, axs=axs, label='Two-infall', marker_labels=True,
                      style_kw={'zorder': 4, 'linestyle': '-.'})
    axs[0].set_xlim((-2.1, 0.6))
    axs[0].set_ylim((-0.2, 0.52))
    axs[0].legend(frameon=False, loc='lower left', handlelength=1.2)
    plt.savefig(paths.figures / 'onezone_sfh.pdf', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()
