r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

from ..._globals import END_TIME
from .utils import exponential
from .normalize import normalize_ifrmode
from .gradient import gradient
import math as m
import os

_FIRST_TIMESCALE = 0.1
_SECOND_TIMESCALE = 4.

class twoinfall:

	def __init__(self, radius, dt = 0.01, dr = 0.1):
		self.first.timescale = 1
		self.second.timescale = 4
		# self.first.timescale = 0.5 + radius / 15.5
		# self.second.timescale = 2 + 4 * (radius / 15.5)
		# self.ratio = self.amp_ratio(radius)
		self.ratio = twoinfall_ampratio(self, gradient, radius,
			onset = self.onset, dt = dt, dr = dr)
		prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
			dr = dr)
		self.first.norm *= prefactor
		self.second.norm *= prefactor