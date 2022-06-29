"""
This file defines the exponential delay-time distribution (DTD) of Type Ia
supernovae.
"""

import math as m

class exponential:
    """
    A normalized, exponentially-declining Type Ia supernova delay time
    distribution.

    Attributes
    ----------
    timescale : float [default: 1.5]
        The exponential timescale in units of Gyr.
    coeff : float [default: 1e-9]
        The post-normalization coefficient. The default is 1e-9 to
        convert between the timescale (in Gyr) and the rate (in yr^-1).

    """
    def __init__(self, timescale=1.5, coeff=1e-9):
        self.timescale = timescale
        self.coeff = coeff

    def __call__(self, time):
        norm = (self.coeff / self.timescale)
        return norm * m.exp(-time / self.timescale)
