"""
This script plots a comparison between one-zone two-infall models with and
without outflows.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import models, dtds
from migration.src._globals import END_TIME
from migration.src.simulations.yields import twoinfall

DT = 0.01
RADII = [4, 6, 8, 10, 12, 14] # kpc

def main(overwrite=True):
    
    output_dir = paths.data / 'onezone' / 'twoinfall'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    fig, axs = plt.subplots(2, 2, figsize=(7, 7), sharex=True)
    
    for galr in RADII:
        area = np.pi * ((galr + 0.1)**2 - galr**2)
        output = str(output_dir / (f'spitoni21{galr:02d}'))
        hist = vice.history(output)
        axs[0,0].plot(hist['time'], np.array(hist['ifr']) / area, label=galr)
        # ifr = models.twoinfall(galr)
        # axs[0,0].plot(simtime, [ifr(t) for t in simtime], c='k', ls='--')
        axs[0,1].plot(hist['time'], np.array(hist['sfr']) / area)
        axs[1,0].plot(hist['time'], np.array(hist['mgas']) / area)
        tau_star = [hist['mgas'][i+1] / hist['sfr'][i+1] * 1e-9 for i in range(
                    len(hist['time']) - 1)]
        axs[1,1].plot(hist['time'][1:], tau_star)
        # tau_star = models.twoinfall_tau_star(area, galr)
        # axs[1,1].plot(hist['time'][1:], [tau_star(hist['time'][i+1], hist['mgas'][i+1]) for i in range(len(hist['time'])-1)], c='k', ls='--')
        
    axs[0,0].set_title('Infall surface density')
    axs[0,1].set_title('Star formation surface density')
    axs[1,0].set_title('Gas mass surface density')
    axs[1,1].set_title(r'$\tau_*$')
    axs[1,1].set_yscale('log')
    axs[1,0].set_xlabel('Time [Gyr]')
    
    axs[0,0].legend(frameon=False)
    fig.savefig(paths.figures / 'onezone_spitoni21_sfh.png', dpi=300)
    

if __name__ == '__main__':
    main()
