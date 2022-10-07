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
from .utils import polyfit


class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        spitoni_params = np.genfromtxt('%s/spitoni_twoinfall.dat' % (
            os.path.abspath(os.path.dirname(__file__))))
        super().__init__(onset = 4, ratio = 0.2) # dummy values
        self.first.timescale = 0.1#polyfit(radius, spitoni_params, 3)
        self.second.timescale = 4#polyfit(radius, spitoni_params, 5)
        self.onset = 4#polyfit(radius, spitoni_params, 9)
        # self.onset = 2 + radius / 4
        # thin_to_thick = polyfit(radius, spitoni_params, 7)
        # self.first.timescale = self.timescale1(radius)
        # self.second.timescale = self.timescale2(radius)
        # self.onset = self.tmax(radius)
        # thin_to_thick = self.thin_to_thick(radius)
        # thin_to_thick = m.exp(radius * (1 / 2 - 1 / 2.5)) / 0.27
        # self.ratio = thin_to_thick * self.timescale_ratio()
        # print(self.ratio)
        self.ratio = 0.1
        # print(self.ratio)
        area = m.pi * ((radius + dr/2) ** 2 - (radius - dr/2) ** 2)
        self.first.norm /= area
        self.second.norm /= area
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr, outflows = True, which_tau_star = 'johnson21')
        # Spitoni's parameters produce an infall rate in terms of mass, but
        # the multizone simulations take surface mass density
        # prefactor /= m.pi * ((radius + dr/2) ** 2 - (radius - dr/2) ** 2)
        # Also, Spitoni's parameters already incorporate the radial gas density
        # gradient
        # prefactor /= gradient(radius)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
    
    def timescale_ratio(self):
        timescale_factor = self.first.timescale / self.second.timescale
        timescale_factor *= (1 - m.exp(-END_TIME / self.first.timescale))
        timescale_factor /= (1 - m.exp(-(END_TIME - self.onset) / 
                                       self.second.timescale))
        return timescale_factor
    
    def timescale1(self, radius):
        if radius <= 8:
            return 0.1
        else:
            return 0.1 * (radius - 7)
        
    def timescale2(self, radius):
        if radius <= 4:
            return 4
        else:
            return 2 * (radius - 6)
        
    def thin_to_thick(self, radius):
        return 3/4 * radius
    
    def tmax(self, radius):
        return -1/4 * radius + 6
