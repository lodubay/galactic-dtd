r"""
This file declares the time-dependence of the infall rate and star formation
efficiency timescale at a given radius in a star formation model which combines
the early history from Conroy et al. (2022) with a late burst from
Johnson et al. (2021).
"""

from ..._globals import END_TIME
from .utils import modified_exponential, gaussian
from .insideout import _TAU_RISE_, insideout
from .normalize import normalize
from .gradient import gradient
import math as m
import os

_BURST_TIME_ = END_TIME - 2 # Gyr

class hybridburst:
	pass
