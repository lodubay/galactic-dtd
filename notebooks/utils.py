# -*- coding: utf-8 -*-
"""
This file contains several basic functional forms which are useful for the
construction of star formation rate, infall rate, and delay time distribution
models.
"""

import math as m

class NormalizedExponential:
    """
    A generic normalized exponentially-declining function of time.

    Attributes
    ----------
    timescale : float
        The exponential timescale in units of time.
    coeff : float
        The post-normalization coefficient.

    """
    def __init__(self, timescale=1, coeff=1):
        self.timescale = timescale
        self.coeff = coeff

    def __call__(self, time):
        return (self.coeff / self.timescale) * m.exp(-time / self.timescale)


class Gaussian:
    """
    A generic normalized Gaussian function of time.

    Attributes
    ----------
    center : float
        The location of the peak of the Gaussian function.
    stdev : float
        The standard deviation of the Gaussian function.
    norm : float
        The normalization of the Gaussian function.

    """
    def __init__(self, center=1, stdev=1, coeff=1, normalize=True):
        """
        Initialize the Gaussian.

        Parameters
        ----------
        center : float [default: 1]
            The location of the peak of the Gaussian function.
        stdev : float [default: 1]
            The standard deviation of the Gaussian function.
        coeff : float [default: 1]
            The post-normalization coefficient.

        """
        self.center = center
        self.stdev = stdev
        self.coeff = coeff
        if normalize:
            self.norm = 1 / (stdev * m.sqrt(2 * m.pi))
        else:
            self.norm = 1

    def __call__(self, time):
        C = self.coeff * self.norm
        return C * m.exp(-(time-self.center)**2 / (2*self.stdev**2))