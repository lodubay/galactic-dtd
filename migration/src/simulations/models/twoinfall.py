r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode
from .gradient import gradient
import math as m
import os

FIRST_TIMESCALE = 0.1 # Gyr
SECOND_TIMESCALE = 4. # Gyr
SECOND_ONSET = 4. # Gyr

THIN_DISK_SCALE_RADIUS = 2.5 # kpc
THICK_DISK_SCALE_RADIUS = 2.0 # kpc
THICK_TO_THIN_RATIO = 0.27 # at r = 0

class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1):
        super().__init__(onset=SECOND_ONSET, ratio=1) # dummy ratio value
        # Calculate the amplitude ratio of infalls
        thin_to_thick = (1 / THICK_TO_THIN_RATIO) * m.exp(radius * (
            1 / THICK_DISK_SCALE_RADIUS - 1 / THIN_DISK_SCALE_RADIUS))
        timescale_ratio = FIRST_TIMESCALE / SECOND_TIMESCALE
        timescale_ratio *= (1 - m.exp(-END_TIME / FIRST_TIMESCALE))
        timescale_ratio /= (1 - m.exp(-(END_TIME - SECOND_ONSET) / SECOND_TIMESCALE))
        self.first.timescale = FIRST_TIMESCALE
        self.second.timescale = SECOND_TIMESCALE
        self.ratio = thin_to_thick * timescale_ratio
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
