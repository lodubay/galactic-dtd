"""
This file defines the bimodal delay-time distribution (DTD) of Type Ia
supernovae.
"""

from .utils import gaussian, exponential
from ..._globals import END_TIME

class bimodal:
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
                 tmin=0.04, tsplit=0.1, tmax=END_TIME):
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
        self.prompt = gaussian(center=center, stdev=stdev, coeff=0.5)
        self.tardy = exponential(timescale=timescale, coeff=0.5)
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
