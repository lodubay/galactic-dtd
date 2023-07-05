"""
This script plots one-zone tracks for the four star formation history (SFH)
models in the paper.
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src import models, dtds
# from multizone.src.yields import C22
from track_and_mdf import setup_axes, plot_vice_onezone
from _globals import END_TIME, DT

def main():
    
    output_dir = paths.data / 'onezone' / 'sfh'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    rgal = 8.
    dr = 0.1
    area = m.pi * ((rgal + dr)**2 - rgal**2)
    
    kwargs = dict(
        mode='sfr',
        elements=('fe', 'o'),
        Mg0=0,
        dt=DT,
        recycling='continuous',
        RIa=dtds.powerlaw(slope=-1.1),
        delay=0.04,
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
                         tau_star=models.earlyburst_tau_star(area))
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
    plot_vice_onezone('../data/onezone/conroy22/conroy22')
    axs[0].set_xlim((-3, 0.5))
    axs[0].set_ylim((-0.1, 0.65))
    plt.savefig(paths.figures / 'onezone_conroy22.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()
