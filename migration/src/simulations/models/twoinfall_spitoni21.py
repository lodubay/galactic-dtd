r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

import numpy as np
from numpy.polynomial.polynomial import Polynomial
from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode, twoinfall_ampratio
from .gradient import gradient
import math as m
import os


class twoinfall_spitoni21(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        spitoni_params = np.genfromtxt('%s/spitoni_twoinfall.dat' % (
            os.path.abspath(os.path.dirname(__file__))))
        super().__init__(onset = 4, ratio = 0.2) # dummy values
        self.first.timescale = self.polyfit(radius, spitoni_params, 3)
        self.second.timescale = self.polyfit(radius, spitoni_params, 5)
        self.onset = self.polyfit(radius, spitoni_params, 9)
        thin_to_thick = self.polyfit(radius, spitoni_params, 7)
        self.ratio = thin_to_thick * self.timescale_ratio()
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr, outflows = False, which_tau_star = 'spitoni21')
        self.first.norm *= prefactor
        self.second.norm *= prefactor
        
    def polyfit(self, radius, params, col):
        fit = Polynomial.fit(params[:,0], params[:,col], deg=2, 
                             w=1/params[:,col+1])
        return fit(radius)
    
    def timescale_ratio(self):
        timescale_factor = self.first.timescale / self.second.timescale
        timescale_factor *= (1 - m.exp(-END_TIME / self.first.timescale))
        timescale_factor /= (1 - m.exp(-(END_TIME - self.onset) / 
                                       self.second.timescale))
        return timescale_factor
