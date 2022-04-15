"""
This file contains classes and functions to specify one-zone model infall rates.
"""

import math as m
from utils import Gaussian
    
class BrokenExponential:
    """
    An infall rate that starts constant and transitions to exponential.
    """
    def __init__(self, tsplit=5, tmax=13.2, timescale=15, coeff=1e-9):
        self.tsplit = tsplit
        self.plateau = coeff
        self.timescale = timescale
        self.exponential = Exponential(timescale=timescale, tmax=tmax, coeff=coeff)
        self.norm = 1
        self.norm = coeff * self.normalize(tmax=tmax)
        
    def __call__(self, time):
        if time < self.tsplit:
            ifr = self.plateau
        else:
            ifr = (self.plateau / m.exp(-self.tsplit/self.timescale)) * m.exp(-time/self.timescale)
        return self.norm * ifr
    
    def normalize(self, tmin=0, tmax=13.2, dt=1e-3):
        integral = 0
        time = tmin
        while time < tmax:
            integral += self.__call__(time) * dt
            time += dt
        return 1 / integral


class Constant:
    """
    A constant normalized infall rate as a function of time.
    """
    def __init__(self, tmax=13.2, coeff=1e-9):
        self.tmax = tmax
        self.coeff = coeff
        
    def __call__(self, time):
        return self.coeff / self.tmax


class Exponential:
    """
    A normalized exponentially declining infall rate.
    """
    def __init__(self, timescale=15, coeff=1e-9, tmax=13.2):
        self.timescale = timescale
        self.coeff = coeff
        self.tmax = tmax

    def __call__(self, time):
        # Normalize the exponential taking the end-point into account
        norm = 1 / ((1 - m.exp(-self.tmax/self.timescale)) * self.timescale)
        return self.coeff * norm * m.exp(-time / self.timescale)


class InsideOut:
    """
    A modification to the exponentially declining IFR with an initial
    exponential rise.
    """
    def __init__(self, rise_timescale=2, ifr_timescale=15):
        self.rise = Exponential(timescale=rise_timescale, coeff=1)
        self.decline = Exponential(timescale=ifr_timescale)
        # Re-normalize to include rise exponential
        self.norm = (1 - 1 / (rise_timescale + ifr_timescale))**-1

    def __call__(self, time):
        return self.norm * (1 - self.rise(time)) * self.decline(time)


class LateBurst(InsideOut):
    """
    A modification to the Inside-Out IFR with the addition of a late,
    Gaussian infall rate burst.
    """
    def __init__(self, burst_strength=1.5, burst_time=10.5, burst_stdev=1,
                 **kwargs):
        self.burst = Gaussian(center=burst_time, stdev=burst_stdev, 
                              coeff=burst_strength)
        super().__init__(**kwargs)
        self.norm *= 1e-9 * self.normalize()

    def __call__(self, time):
        return (1 + self.burst(time)) * super().__call__(time)
    
    def normalize(self, tmin=0, tmax=13.2, dt=1e-3):
        integral = 0
        time = tmin
        while time < tmax:
            integral += self.__call__(time) * dt
            time += dt
        return 1 / integral