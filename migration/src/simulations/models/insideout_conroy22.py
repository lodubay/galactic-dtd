"""
This file declares a modified version of the inside-out star formation history
for use with the conroy22 SFE timescale model.
"""

from .insideout import insideout
from .normalize import normalize_conroy22
from .gradient import gradient


class insideout_conroy22(insideout):
	def __init__(self, radius, dt = 0.01, dr = 0.01):
		super().__init__(radius, dt = dt, dr = dr)
		self.norm = normalize_conroy22(self, gradient, radius, dt = dt, dr = dr)
