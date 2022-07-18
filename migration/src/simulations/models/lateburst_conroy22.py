"""
This file declares a modified version of the late-burst star formation history
for use with the conroy22 SFE timescale model.
"""

from .lateburst import lateburst
from .normalize import normalize_conroy22
from .gradient import gradient


class lateburst_conroy22(lateburst):
    def __init__(self, radius, dt = 0.01, dr = 0.01):
        super().__init__(radius, dt = dt, dr = dr)
        self._prefactor = normalize_conroy22(self, gradient, radius, dt = dt,
                                             dr = dr)
