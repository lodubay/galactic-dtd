"""
This script plots a comparison between one-zone two-infall models with and
without outflows.
"""

import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from track_and_mdf import setup_axes, plot_vice_onezone
from multizone.src.simulations import models, dtds
from _globals import END_TIME
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
from vice.toolkit import J21_sf_law

DT = 0.01
RADII = [4, 6, 8, 10, 12] # kpc

def main(overwrite=True):
    
    output_dir = paths.data / 'onezone' / 'insideout'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        
    kwargs = dict(
        mode='sfr',
        elements=('fe', 'o'),
        Mg0=0,
        dt=DT,
        recycling='continuous',
        RIa=dtds.powerlaw(slope=-1.1),
        delay=0.04,
        schmidt=False
    )
    simtime = np.arange(0, END_TIME+DT, DT)
    
    fig, axs = setup_axes()
    
    for galr in RADII:
        area = np.pi * ((galr + 0.1)**2 - galr**2)
        # print(J21_sf_law(area)(0, 0))
        output = str(output_dir / (f'insideout{galr:02d}'))
        sz = vice.singlezone(name=output,
                             tau_star=2,
                             func=models.insideout(galr, dt=DT),
                             eta=vice.milkyway.default_mass_loading(galr),
                             **kwargs)
        sz.run(simtime, overwrite=overwrite)
        # hist = vice.history(output)
        # print(hist)
        plot_vice_onezone(output, fig=fig, axs=axs, label='%s kpc' % galr,
                          marker_labels=(galr==4.))
        
    axs[0].legend(loc='lower left', frameon=False)
    axs[0].set_xlim((-2.5, 0.9))
    axs[0].set_ylim((-0.24, 0.54))
    fig.savefig(paths.figures / 'onezone_insideout.png', dpi=300)
    

if __name__ == '__main__':
    main()
