"""
This script plots one-zone tracks for the four star formation history (SFH)
models in the paper.
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
from multizone.src import models, dtds
# from multizone.src.yields import C22
from track_and_mdf import setup_axes, plot_vice_onezone
from _globals import END_TIME, DT, MIN_RIA_DELAY

def main():
    
    output_dir = paths.data / 'onezone' / 'sfh'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    rgal = 8.
    dr = 0.1
    area = m.pi * ((rgal + dr)**2 - rgal**2)
    
    kwargs = dict(
        elements=('fe', 'o'),
        Mg0=0,
        dt=DT,
        recycling='continuous',
        RIa=dtds.powerlaw(slope=-1.1),
        delay=MIN_RIA_DELAY,
        schmidt=False,
        eta=vice.milkyway.default_mass_loading(rgal)
    )
    simtime = np.arange(0, END_TIME+DT, DT)
    
    # Inside-out SFH
    sz = vice.singlezone(name=str(output_dir / 'insideout'),
                         func=models.insideout(rgal, dr=dr), 
                         mode='sfr',
                         tau_star=vice.toolkit.J21_sf_law(area, mode='sfr'),
                         **kwargs)
    sz.run(simtime, overwrite=True)
    
    # Late-burst SFH
    sz = vice.singlezone(name=str(output_dir / 'lateburst'),
                         func=models.lateburst(rgal, dr=dr), 
                         mode='sfr',
                         tau_star=vice.toolkit.J21_sf_law(area, mode='sfr'),
                         **kwargs)
    sz.run(simtime, overwrite=True)
    
    # Early-burst SFH
    sz = vice.singlezone(name=str(output_dir / 'earlyburst'),
                         func=models.earlyburst_ifr(rgal, dr=dr), 
                         mode='ifr',
                         tau_star=models.earlyburst_tau_star(area),
                         **kwargs)
    sz.run(simtime, overwrite=True)
    
    # Two-infall SFH
    sz = vice.singlezone(name=str(output_dir / 'twoinfall'),
                         func=models.twoinfall(rgal, dr=dr), 
                         mode='ifr',
                         tau_star=vice.toolkit.J21_sf_law(area, mode='ifr'),
                         **kwargs)
    sz.run(simtime, overwrite=True)
    
    # Plot
    fig, axs = setup_axes()
    plot_vice_onezone(str(output_dir / 'insideout'), fig=fig, axs=axs, label='Inside-out')
    plot_vice_onezone(str(output_dir / 'lateburst'), fig=fig, axs=axs, label='Late-burst')
    plot_vice_onezone(str(output_dir / 'earlyburst'), fig=fig, axs=axs, label='Early-burst')
    plot_vice_onezone(str(output_dir / 'twoinfall'), fig=fig, axs=axs, label='Two-infall', marker_labels=True)
    axs[0].set_xlim((-1.6, 0.6))
    axs[0].set_ylim((-0.2, 0.52))
    plt.savefig(paths.figures / 'onezone_sfh.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()
