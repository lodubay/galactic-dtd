r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode, twoinfall_ampratio
from .gradient import gradient
import math as m
# import os

FIRST_TIMESCALE = 1. # Gyr
SECOND_TIMESCALE = 4. # Gyr
SECOND_ONSET = 4. # Gyr, at R = 8 kpc

# THIN_DISK_SCALE_RADIUS = 2.5 # kpc
# THICK_DISK_SCALE_RADIUS = 2.0 # kpc
# THIN_TO_THICK_RATIO = 5.6 # at R = 0 kpc

class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        super().__init__() # dummy initial parameters
        self.onset = SECOND_ONSET #+ (radius - 8)/4
        # self.onset = SECOND_ONSET
        # Calculate the amplitude ratio of infalls
        # thin_to_thick = THIN_TO_THICK_RATIO * m.exp((radius - 8) * (
        #     1 / THICK_DISK_SCALE_RADIUS - 1 / THIN_DISK_SCALE_RADIUS))
        # timescale_ratio = FIRST_TIMESCALE / SECOND_TIMESCALE
        # timescale_ratio *= (1 - m.exp(-END_TIME / FIRST_TIMESCALE))
        # timescale_ratio /= (1 - m.exp(-(END_TIME - self.onset) / SECOND_TIMESCALE))
        self.first.timescale = FIRST_TIMESCALE 
        # self.first.timescale = max(0.1, 0.1 + (radius - 8) * 0.1)
        self.second.timescale = SECOND_TIMESCALE 
        # self.second.timescale = max(4, 4 + (radius - 8) * 1.5)
        # self.ratio = thin_to_thick * timescale_ratio
        self.ratio = twoinfall_ampratio(self, radius, onset=self.onset, 
                                        dr = dr, dt = dt)
        # self.ratio = 1 / 3.5
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
