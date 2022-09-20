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

def main(evolutions, RIa='powerlaw_slope11', migration='diffusion'):
    fig, ax = plt.subplots()
    
    for evolution in evolutions:
        output_name = 'migration/%s/%s/%s' % (migration, evolution, RIa)
        multiout = vice.multioutput(str(paths.data / output_name))
        # print(multiout.zones['zone1'].history["mstar"][-1])
        for zone in multiout.zones:
            print(zone)
        # stars = multioutput_to_pandas(output_name)
        
        radii = np.arange(0, 20+ZONE_WIDTH, ZONE_WIDTH)
        sigma_star = [zone.history["mstar"][-1] for zone in multiout.zones]
        print(sigma_star)
        # sigma_star = np.zeros(radii.shape)
        # for i, radius in enumerate(radii):
        #     area = np.pi * ((radius + ZONE_WIDTH)**2 - radius**2)
        #     subset = filter_multioutput_stars(stars, 
        #                                       galr_lim=(radius, radius + ZONE_WIDTH),
        #                                       zone_width=ZONE_WIDTH)
        #     sigma_star[i] = subset['mass'].sum() / area
        
        ax.plot(radii, sigma_star, label=evolution)
    
    ax.set_xlabel('Galactocentric radius [kpc]')
    ax.set_ylabel(r'$\Sigma_*$')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(paths.figures / 'stellar_mass_surface_density.png', dpi=300)

if __name__ == "__main__":
    evolutions = sys.argv[1:]
    main(evolutions)
