"""
This file contains the Type Ia supernova delay time distributions and related
functions.
"""

import math as m
from utils import Gaussian

class Bimodal:
    """
    The bimodal delay-time distribution of SNe Ia. This assumes ~50% of SNe Ia
    belong to a prompt (<0.1 Gyr) component with the form of a narrow Gaussian,
    and the remaining ~50% form an exponential DTD.

    Attributes
    ----------
    prompt : <function>
        Prompt component as a function of time.
    tardy : <function>
        Tardy component as a function of time.
    norm : float
        Normalization coefficient scaled so the total integral is unity.
    """

    def __init__(self, center=0.05, stdev=0.01, timescale=3,
                 tmin=0.04, tsplit=0.1, tmax=13.2):
        """
        Initialize the bimodal model.

        Parameters
        ----------
        center : float [default: 0.05]
            Center of the prompt Gaussian component in Gyr.
        stdev : float [default: 0.01]
            Standard deviation of the prompt Gaussian component in Gyr.
        timescale : float [default: 3]
            Exponential timescale of the tardy component in Gyr.
        tmin : float [default: 0.04]
            Minimum delay time in Gyr for integration purposes.
        tsplit : float [default: 0.1]
            Time in Gyr to switch from prompt to tardy component.
        tmax : float [default: 13.2]
            Maximum delay time in Gyr for integration purposes.

        """
        self.prompt = Gaussian(center=center, stdev=stdev, coeff=0.5)
        self.tardy = Exponential(timescale=timescale, coeff=0.5)
        self.norm = 1
        # Normalize over full time range
        self.norm *= 1e-9 * self.normalize(tmin, tmax)

    def __call__(self, time):
        """
        Calculate the normalized SN Ia rate at the given time.

        Parameters
        ----------
        time : float
            Time in Gyr since the starburst.

        Returns
        -------
        RIa : float
            Normalized SN Ia rate per solar mass per year.

        """
        return self.norm * (self.prompt(time) + self.tardy(time))

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


class BrokenPowerLaw:
    """
    A two-part broken power-law delay-time distribution of SNe Ia. The default
    setting is a flat distribution (slope of 0) before the time of separation
    and the standard -1.1 slope after.

    Attributes
    ----------
    tsplit : float
        Time in Gyr separating the two power-law components.
    plaw1 : <function>
        The first power-law function, called when t < tsplit.
    plaw2 : <function>
        The second power-law function, called when t >= tsplit.

    """
    def __init__(self, tsplit=0.2, slope1=0, slope2=-1.1, tmin=0.04,
                 tmax=13.2):
        """
        Initialize the broken power-law.

        Parameters
        ----------
        tsplit : float [default: 0.2]
            Time in Gyr separating the two components.
        slope1 : float [default: 0]
            The slope of the first power-law.
        slope2 : float [default: -1.1]
            The slope of the second power-law.
        tmin : float [default: 0.04]
            The minimum delay time and lower limit of integration.
        tmax : float [default: 13.2]
            The maximum simulation time and upper limimt of integration.

        """
        self.tsplit = tsplit
        # Initialize both power-law functions
        plaw1 = PowerLaw(slope=slope1, tmin=tmin, tmax=tsplit, coeff=1e-9)
        plaw2 = PowerLaw(slope=slope2, tmin=tsplit, tmax=tmax, coeff=1e-9)
        # Calculate new normalization coefficients
        norm1 = (plaw1.norm**-1 + tsplit**(slope1-slope2) * plaw2.norm**-1)**-1
        plaw1.norm = norm1
        norm2 = tsplit**(slope1-slope2) * norm1
        plaw2.norm = norm2
        self.plaw1 = plaw1
        self.plaw2 = plaw2

    def __call__(self, time):
        """
        Calculate the normalized SN Ia rate at a given time.

        Parameters
        ----------
        time : float
            Time since starburst in Gyr.

        Returns
        -------
        RIa : float
            Normalized SN Ia rate per solar mass per year.

        """
        if time < self.tsplit:
            RIa = self.plaw1(time)
        else:
            RIa = self.plaw2(time)
        return RIa


class Exponential:
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


class PowerLaw:
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
    def __init__(self, slope=-1.1, coeff=1e-9, tmin=0.04, tmax=13.2):
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


