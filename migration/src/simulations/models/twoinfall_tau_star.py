r"""
This file implements a version of the star formation efficiency model from 
Spitoni et al. (2021) to be used in VICE multizone simulations with the
two-infall star formation history.
"""

import os
import numpy as np
from vice.toolkit import J21_sf_law
from .utils import polyfit

class twoinfall_tau_star(J21_sf_law):
    r"""
    An implementation of the Spitoni et al. (2021) SFE timescale model.
    
    This class maintains the dependence of $\tau_*$ on $M_{\rm{gas}}$ as a
    three-component power-law which is implemented in `vice.toolkit.J21_sf_law`.
    However, it strips out the time dependence (based on the SFE timescale
    of molecular gas as a function of redshift per Tacconi et al. 2018) and
    replaces it with the model by Conroy et al. (2021), also introducing a
    radial dependence.
    """
    def __init__(self, area, radius, mode = "ifr", **kwargs):
        spitoni_params = np.genfromtxt('%s/spitoni_twoinfall.dat' % (
            os.path.abspath(os.path.dirname(__file__))))
        self.radius = radius
        self.onset = 4#polyfit(radius, spitoni_params, 9)
        super().__init__(area, mode = mode, **kwargs)
        
    def __call__(self, time, mgas):
        mgas_dependence = super().__call__(time, mgas) / self.molecular(time)
        return mgas_dependence * self.time_dependence(time)
        # return self.time_dependence(time)
        # return mgas_dependence

    def time_dependence(self, time):
        r"""
        Implementation of the star formation efficiency timescale prescription
        in Table 2 of Spitoni et al. (2021).

        Parameters
        ----------
        time : float
            Time after the starburst in Gyr.

        Returns
        -------
        float
            Star formation efficiency timescale $\tau_*$ in Gyr.

        """
        if time < self.onset:
            return 1 / self.sfe1
        else:
            return 1 / self.sfe2
      
    @property
    def sfe1(self):
        """
        float
            The SFE of the first (high-alpha) infall in Gyr^-1.
        """
        if self.radius < 8.:
            return 4 - (self.radius / 4)
        else:
            return 2
        
    @property
    def sfe2(self):
        """
        float
            The SFE of the second (low-alpha) infall in Gyr^-1.
        """
        return max(2 - (self.radius / 8), 1e-6)
