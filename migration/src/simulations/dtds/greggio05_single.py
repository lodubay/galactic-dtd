"""
This file contains the greggio05_single class, which implements the single-
degenerate Type Ia DTD from Greggio 2005.
"""

from numbers import Number
import numpy as np
import vice
from ..._globals import END_TIME
from .utils import mlr_wrapper, minimum_wd_mass

class greggio05_single:
    """
    An implementation of the analytic single-degenerate DTD from Greggio 2005.

    Attributes
    ----------
    m2_slope : float
        The power-law slope of the distribution of secondary masses
    q_slope : float
        The power-law slope of the distribution of mass ratios
    efficiency : float
        The mass transfer efficiency
    mlr : str
        Which mass-lifetime relation available in VICE to assume
    imf : str
        Which initial mass function to assume

    Methods
    -------
    secondary_mass_distribution(self, m2, m1_max=8., dm=0.01)
        The distribution of secondary masses in SD progenitor systems
    normalize(self, tmin=0.04, tmax=END_TIME, dt=1e-3)
        Normalize the DTD to unity
    """
    def __init__(self, m2_slope=-1.44, q_slope=1, efficiency=1.,
                 mlr='larson1974', imf='kroupa'):
        """
        Parameters
        ----------
        m2_slope : float, optional
            Slope of the power-law derivative of secondary mass. The default is
            -1.44
        q_slope : float, optional
            The power-law slope of the mass ratio distribution q = m2/m1.
            The default is 1.
        efficiency : float, optional
            Mass-transfer efficiency represented by epsilon in Equation 17.
            The default is 1.
        mlr : str, optional
            Which mass-lifetime relation to use. The default is 'larson1974'
        imf : str, optional
            Which IMF to use. The default is 'kroupa'
        """
        self.m2_slope = m2_slope
        self.q_slope = q_slope
        self.efficiency = efficiency
        self.mlr = mlr
        self.imf = imf
        self.norm = 1
        self.norm = self.normalize()

    def __call__(self, time):
        """
        Compute the single-degenerate DTD at the given time

        Parameters
        ----------
        time : float
            Time since star formation event in Gyr

        Returns
        -------
        float
            The value of the single-degenerate DTD
        """
        secondary_mass = mlr_wrapper(time, which='age', model=self.mlr)
        f_m2 = self.secondary_mass_distribution(secondary_mass)
        return self.norm * f_m2 * time ** self.m2_slope

    @property
    def m2_slope(self):
        """
        Type: float
            The power-law slope of the derivative of the secondary mass
            function
        """
        return self._m2_slope

    @m2_slope.setter
    def m2_slope(self, value):
        if isinstance(value, Number):
            self._m2_slope = value
        else:
            raise TypeError('Parameter "m2_slope" must be numeric. Got: %s' \
                            % type(value))

    @property
    def q_slope(self):
        """
        Type: float
            The power-law slope of the mass ratio distribution (q = m2/m1)
        """
        return self._q_slope

    @q_slope.setter
    def q_slope(self, value):
        if isinstance(value, Number):
            self._q_slope = value
        else:
            raise TypeError('Parameter "q_slope" must be numeric. Got: %s' \
                            % type(value))

    @property
    def efficiency(self):
        """
        Type: float
            The mass-transfer efficiency represented by epsilon in Equation 17
            from Greggio 2005.
        """
        return self._efficiency

    @efficiency.setter
    def efficiency(self, value):
        if isinstance(value, Number):
            if value >= 0 and value <= 1:
                self._efficiency = float(value)
            else:
                raise ValueError('Efficiency must be between 0 and 1.')
        else:
            raise TypeError('Parameter "efficiency" must be numeric. Got: %s' \
                            % type(value))

    @property
    def mlr(self):
        r"""
        Type: str
            Which formulation of the mass-lifetime relation (MLR) to assume.
            Allowed values:

                - 'larson1974'
                - 'mm1989'
                - 'pm1993'
                - 'ka1997'
                - 'hpt2000'
                - 'vincenzo2016'
                - 'powerlaw'
        """
        return self._mlr

    @mlr.setter
    def mlr(self, value):
        if isinstance(value, str):
            try:
                mlr_wrapper(1, model=value)
                self._mlr = value
            except:
                raise ValueError('Parameter "mlr" not in acceptable list.')
        else:
            raise TypeError('Parameter "mlr" must be a string. Got: %s' \
                            % type(value))

    @property
    def imf(self):
        r"""
        Type: str
            Which formulation of the initial mass function (IMF) to assume.
            Allowed values:

                - 'kroupa'
                - 'salpeter'
        """
        return self._imf

    @imf.setter
    def imf(self, value):
        if isinstance(value, str):
            if value.lower() in ['kroupa', 'salpeter']:
                self._imf = value
            else:
                raise ValueError('IMF must be either "kroupa" or "salpeter".')
        else:
            raise TypeError('Parameter "imf" must be a string. Got: %s' \
                            % type(value))

    def secondary_mass_distribution(self, m2, m1_max=8., dm=0.01):
        """
        The distribution function of secondary masses in SD systems.

        Taken from Equation 16 in Greggio 2005. This assumes the mass of the
        primary follows a Kroupa or Salpeter distribution, and the mass ratio
        q = m2/m1 follows a power-law with a slope of 'q_slope'.

        Parameters
        ----------
        m2 : float
            Mass of the secondary in solar masses
        m1_max : float
            Maximum primary mass in solar masses. The default is 8.
        dm : float
            Mass integration step in solar masses. The default is 0.01

        Returns
        -------
        float
            The relative frequency of the given secondary mass
        """
        # Lower limit of integration over primary mass
        m_wd_min = minimum_wd_mass(m2, efficiency=self.efficiency)
        m1_min = max((m2, 2, 2 + 10 * (m_wd_min - 0.6)))
        primary_masses = np.arange(m1_min, m1_max + dm, dm)
        imf_func = {
            'kroupa': vice.imf.kroupa,
            'salpeter': vice.imf.salpeter
        }[self.imf]
        integral = np.sum([m2 ** self.q_slope * m1 ** (-(self.q_slope + 1)) * \
                           imf_func(m1) * dm \
                           for m1 in primary_masses])
        return integral

    def normalize(self, tmin=0.04, tmax=END_TIME, dt=0.01):
        """
        Normalize the area under the DTD to 1.

        Parameters
        ----------
        tmin : float [default: 0.04]
            Minimum integration time in Gyr
        tmax : float [default: 13.2]
            Maximum integration time in Gyr
        dt : float [default: 1e-3]
            Integration timestep in Gyr

        Returns
        -------
        float
            Normalization factor in Msun^-1 yr^-1
        """
        tarr = np.arange(tmin, tmax+dt, dt)
        integral = np.sum([self.__call__(t) * dt * 1e9 for t in tarr])
        return 1 / integral
