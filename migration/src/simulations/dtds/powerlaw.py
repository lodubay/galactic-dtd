"""
This file defines the power-law delay-time distribution (DTD) of Type Ia
supernovae.
"""

from ..._globals import END_TIME

class powerlaw:
    """
    A normalized power-law delay time distribution.

    Attributes
    ----------
    slope : float [default: -1.1]
        The slope of the power-law.
    coeff : float [default: 1e-9]
        The post-normalization coefficient. The default is 1e-9 to
        convert between the timescale (in Gyr) and the rate (in yr^-1).
    norm : float
        The normalization coefficient, determined by integrating over the
        given range.

    """
    def __init__(self, slope=-1.1, coeff=1e-9, tmin=0.04, tmax=END_TIME):
        """
        Initialize the power-law function.

        Parameters
        ----------
        tmin : float [default: 0.04]
            The lower bound in Gyr of range over which to normalize.
        tmax : float [default: 13.2]
            The upper bound in Gyr of range over which to normalize.

        """
        self.slope = slope
        self.coeff = coeff
        self.norm = self.normalize(tmin, tmax)

    def __call__(self, time):
        return self.coeff * self.norm * (time ** self.slope)

    def normalize(self, tmin, tmax):
        intslope = self.slope + 1 # The slope of the integral
        return intslope / (tmax ** intslope - tmin ** intslope)
