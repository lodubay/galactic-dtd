"""
This script plots a comparison between one-zone two-infall models with and
without outflows.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import vice
from vice.yields.presets import JW20
vice.yields.sneia.settings['fe'] *= 10**0.1
import paths
from track_and_mdf import setup_axes, plot_vice_onezone
sys.path.append(str(paths.root))
from migration.src.simulations import models, dtds
from migration.src._globals import END_TIME

MASS_LOADING = [0, 0.5, 1, 2]
DT = 0.01
GALR = 8.0 # kpc

def main(overwrite=True):
    
    output_dir = paths.data / 'onezone' / 'twoinfall'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        
    kwargs = dict(
        func=models.twoinfall(GALR, dt=DT),
        mode='ifr',
        Mg0=0,
        elements=('fe', 'o'),
        dt=DT,
        recycling='continuous',
        RIa=dtds.powerlaw(slope=-1.1),
        delay=0.04,
        tau_star=2
    )
    simtime = np.arange(0, END_TIME+DT, DT)
    
    fig, axs = setup_axes()
    
    for eta in MASS_LOADING:
        output = str(output_dir / ('eta{:02d}'.format(int(eta * 10))))
        sz = vice.singlezone(name=output, eta=eta, **kwargs)
        sz.run(simtime, overwrite=overwrite)
        plot_vice_onezone(output, fig=fig, axs=axs, label=rf'$\eta={eta:.1f}$')
        
    axs[0].legend(loc='lower left', frameon=False)
    axs[0].set_xlim((-2.5, 1))
    axs[0].set_ylim((-0.2, 0.5))
    fig.savefig(paths.figures / 'onezone_twoinfall.png', dpi=300)    

if __name__ == '__main__':
    main()
