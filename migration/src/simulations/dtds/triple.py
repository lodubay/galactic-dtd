"""
This file defines a delay time distribution which approximates estimates from
triple system evolution by Rajamuthukumar et al. (2022).
"""

import math as m
from ..._globals import END_TIME

LOG_T_MIN = -1.
LOG_T_RISE = -0.3
LOG_T_PEAK = 0.
LOG_PLATEAU = -2.25
LOG_PEAK = -1.
DECLINE_SLOPE = -1.

class triple:
    """
    An approximation of the Rajamuthukumar et al. (2022) triple-system DTD.
    
    This class defines a delay time distribution characterized by an initial
    (low) plateau, a rapid rise between 500 Myr and 1 Gyr, and a t^-1 declining
    power-law after 1 Gyr.
    """
    def __init__(self, tmin=10**LOG_T_MIN, tmax=END_TIME):
        self.trise = 10 ** LOG_T_RISE
        self.tpeak = 10 ** LOG_T_PEAK
        self.norm = 1
        # Normalize over full time range
        self.norm *= 1e-9 * self.normalize(tmin, tmax)
        self._name = 'triple_delay{:03d}'.format(int(tmin * 1000))
    
    def __call__(self, time):
        R_plateau = 10 ** LOG_PLATEAU
        R_peak = 10 ** LOG_PEAK
        if time <= self.trise:
            return self.norm * R_plateau 
        elif time <= self.tpeak:
            timescale = -self.tpeak / m.log(1e-4) # small number to bring the exponential rise as close as possible to the peak
            return self.norm * (R_plateau + (R_peak - R_plateau) * (1 - m.exp(-(time - self.trise) / timescale)))
            # rise_slope = (LOG_PEAK - LOG_PLATEAU) / (LOG_T_PEAK - LOG_T_RISE)
            # return self.norm * 10 ** LOG_PEAK * (time / self.tpeak) ** rise_slope
        else:
            return self.norm * R_peak * (time / self.tpeak) ** DECLINE_SLOPE

    def normalize(self, tmin, tmax, dt=1e-3):
        """
        Calculate the normalization coefficient over the whole distribution.
    
        Parameters
        ----------
        tmin : float
            Lower limit of integration in Gyr.
        tmax : float
            Upper limit of integration in Gyr.
        dt : float [default: 1e-3]
            The numerical integration timestep in Gyr.
    
        Returns
        -------
        norm : float
            Normalization coefficient for which the total integral is unity.
    
        """
        integral = 0
        time = tmin
        while time < tmax:
            integral += self.__call__(time) * dt
            time += dt
        return 1 / integral
    