r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode, twoinfall_ampratio
from .gradient import gradient
import math as m

FIRST_TIMESCALE = 1. # Gyr
SECOND_TIMESCALE = 4. # Gyr
SECOND_ONSET = 4. # Gyr, at R = 8 kpc

class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        super().__init__() # dummy initial parameters
        self.onset = SECOND_ONSET
        self.first.timescale = FIRST_TIMESCALE 
        self.second.timescale = SECOND_TIMESCALE 
        self.ratio = twoinfall_ampratio(self, radius, onset=self.onset, 
                                        dr = dr, dt = dt)
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
