r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

import math as m
import os
import numpy as np
from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode
from .gradient import gradient
from .utils import polyfit, expfit


class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        spitoni_params = np.genfromtxt('%s/spitoni_twoinfall.dat' % (
            os.path.abspath(os.path.dirname(__file__))))
        super().__init__(onset = 4, ratio = 0.2) # dummy values
        self.first.timescale = polyfit(radius, spitoni_params, 3)
        self.second.timescale = polyfit(radius, spitoni_params, 5)
        self.onset = polyfit(radius, spitoni_params, 9)
        thin_to_thick = polyfit(radius, spitoni_params, 7)
        self.ratio = thin_to_thick * self.timescale_ratio()
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr, outflows = False, which_tau_star = 'spitoni21')
        # Spitoni's parameters produce an infall rate in terms of mass, but
        # the multizone simulations take surface mass density
        prefactor /= m.pi * ((radius + dr/2) ** 2 - (radius - dr/2) ** 2)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
    
    def timescale_ratio(self):
        timescale_factor = self.first.timescale / self.second.timescale
        timescale_factor *= (1 - m.exp(-END_TIME / self.first.timescale))
        timescale_factor /= (1 - m.exp(-(END_TIME - self.onset) / 
                                       self.second.timescale))
        return timescale_factor
