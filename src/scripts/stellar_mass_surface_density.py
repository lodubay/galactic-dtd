r"""
Plot the stellar mass surface density ($\Sigma_*$) as a function of radius.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars
from _globals import ZONE_WIDTH

def main(evolutions, RIa='exponential_timescale15', migration='post-process'):
    fig, ax = plt.subplots()
    
    for evolution in evolutions:
        output_name = 'migration/%s/%s/%s' % (migration, evolution, RIa)
        multiout = vice.multioutput(str(paths.data / output_name))
        
        radii = np.arange(0, 15, ZONE_WIDTH)
        areas = np.pi * ((radii + ZONE_WIDTH)**2 - radii**2)
        dt = multiout.zones['zone0'].history['time'][1] - multiout.zones['zone0'].history['time'][0]
        sigma_star = [multiout.zones['zone%s' % i].history['mstar'][-1] / areas[i] for i in range(len(radii))]
        
        # integrated_sfr = np.zeros(radii.shape)
        # for i in range(len(radii)):
        #     zone = multiout.zones['zone%s' % i]
        #     integrated_sfr[i] = sum(zone.history['sfr']) * dt * 1e9
        # integrated_sf_density = integrated_sfr / areas
        
        if evolution in ['insideout', 'expifr_outflows', 'expifr_no_outflows']:
            ls = '-'
        else:
            ls = '--'
        ax.plot(radii, sigma_star, label=evolution, ls=ls)
        # ax.plot(radii, integrated_sf_density, ls='--')
    
    ax.set_xlabel('Galactocentric radius [kpc]')
    ax.set_ylabel(r'$\Sigma_*$ [$M_{\odot}\,\rm{kpc}^{-2}$]')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(paths.figures / 'stellar_mass_surface_density.png', dpi=300)

if __name__ == "__main__":
    evolutions = sys.argv[1:]
    main(evolutions)
