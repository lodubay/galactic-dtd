"""
This script plots a comparison between one-zone two-infall models with and
without outflows.
"""

import math as m
import sys
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from track_and_mdf import setup_axes, plot_vice_onezone
sys.path.append(str(paths.root))
from migration.src.simulations import models, dtds
from migration.src._globals import END_TIME
from migration.src.simulations.yields import twoinfall
from migration.src.simulations.models.gradient import gradient

ZONE_WIDTH = 0.1
DT = 0.01
RADII = [4, 6, 8, 10, 12, 14] # kpc

class twoinfall_gradient(models.twoinfall):
    def __init__(self, radius, dt=DT, dr=ZONE_WIDTH):
        super().__init__(radius, dt=dt, dr=dr)
        area = m.pi * ((radius + dr/2)**2 - (radius - dr/2)**2)
        self.first.norm *= area * gradient(radius)
        self.second.norm *= area * gradient(radius)

def main(overwrite=True):
    
    output_dir = paths.data / 'onezone' / 'twoinfall'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        
    kwargs = dict(
        mode='ifr',
        Mg0=0,
        eta=0,
        elements=('fe', 'o'),
        dt=DT,
        recycling='continuous',
        RIa=dtds.exponential(timescale=1.5),
        delay=0.04,
        schmidt=False
    )
    simtime = np.arange(0, END_TIME+DT, DT)
    
    fig, axs = setup_axes()
    
    for galr in RADII:
        area = np.pi * ((galr + ZONE_WIDTH)**2 - galr**2)
        galr_mean = galr + ZONE_WIDTH / 2
        output = str(output_dir / (f'twoinfall{galr:02d}'))
        sz = vice.singlezone(name=output,
                             tau_star=models.twoinfall_tau_star(area, galr_mean),
                             func=twoinfall_gradient(galr_mean),
                             # eta=vice.milkyway.default_mass_loading(galr),
                             **kwargs)
        sz.run(simtime, overwrite=overwrite)
        plot_vice_onezone(output, fig=fig, axs=axs, label='%s kpc' % galr,
                          marker_labels=(galr==4.))
        
    axs[0].legend(loc='lower left', frameon=False)
    axs[0].set_xlim((-2.5, 0.9))
    axs[0].set_ylim((-0.24, 0.54))
    fig.savefig(paths.figures / 'onezone_twoinfall.png', dpi=300)
    

if __name__ == '__main__':
    main()
