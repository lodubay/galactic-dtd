"""
This file contains a prescription for the DTD of SNe Ia following the double-
degenerate solution of Greggio 2005.
"""

from numbers import Number
import math as m
import numpy as np
from tqdm import tqdm
import vice
from ..._globals import END_TIME

class greggio05_double:
    def __init__(self, scheme='wide', t_nuc_min=0.04, t_nuc_max=1.,
                 t_grav_min=1e-3, beta_sep=0., beta_grav=-0.75, dt=1e-3,
                 mlr='larson1974', nsamples=100, efficiency=1., **kwargs):
        """
        Parameters
        ----------
        scheme : str [default: 'wide']
            Keyword denoting assumptions about the evolution of a DD system
        t_nuc_min : float [default: 0.04]
            Nuclear timescale in Gyr of the most massive possible secondary,
            assumed to be 8 solar masses
        t_nuc_max : float [default: 1]
            Nuclear timescale in Gyr of the least massive possible secondary,
            assumed to be 2 solar masses
        t_grav_min : float [default: 1e-3]
            Minimum gravitational delay in Gyr of the DD progenitor, assumed
            to be independent of the mass of the secondary
        beta_sep : float [default: 0]
            The exponent of the power law distribution adopted for the final
            separations in the 'wide' scheme; a flat distribution of
            of separations corresponds to beta_sep = 0
        beta_grav : float [default: -0.75]
            The exponent of the power law distribution adopted for the
            gravitational delays in the 'close' scheme; a flat distribution
            of final separations corresponds to beta_grav = -0.75
        dt : float [default: 1e-3]
            Integration timestep in Gyr
        mlr : str [default: 'larson1974']
            Which model of mass-lifetime relation to use. Must be one available
            in VICE
        nsamples : int [default: 100]
            Number of log-spaced points at which to pre-compute the DTD
        efficiency : float [default: 0.5]
            Mass-transfer efficiency represented by epsilon in Equation 17
        Other keyword arguments are passed to greggio05_single
        """
        self.scheme = scheme
        self.t_nuc_min = check_if_instance(t_nuc_min, float, 't_nuc_min')
        self.t_nuc_max = check_if_instance(t_nuc_max, float, 't_nuc_max')
        self.t_grav_min = check_if_instance(t_grav_min, float, 't_grav_min')
        self.dt = dt
        # Total minimum and maximum delay times
        self.t_min = t_nuc_min + t_grav_min
        self.t_max = t_nuc_max + self.maximum_gravitational_delay(t_nuc_max)

        if scheme == 'wide':
            self.beta = beta_sep
        else:
            self.beta = beta_grav

        self.efficiency = efficiency
        self.mlr = mlr

        self.single_degenerate_distribution = greggio05_single(
            efficiency=self.efficiency, mlr=self.mlr, **kwargs
        )

        # Pre-calculate the DTD to save on compute time
        print('Computing Greggio 2005 DTD...')
        self.times = np.logspace(np.log10(self.t_min),
                                 np.log10(min((END_TIME, self.t_max))),
                                 num=nsamples, endpoint=True)
        self.dtd = np.zeros(self.times.shape)
        for i, t in enumerate(tqdm(self.times)):
            self.dtd[i] = self.integrate(t)

        # Normalize the DTD
        timesteps = self.times[1:] - self.times[:-1]
        self.dtd /= np.sum(self.dtd[:-1] * timesteps * 1e9)

    def __call__(self, time):
        """
        Linearly interpolate the DTD at the given time.

        Parameters
        ----------
        time : float
            Time after starburst in Gyr

        Returns
        -------
        float
            Interpolated value of the distribution function of DD SNe Ia,
            or NaN if time is outside the bounds of integration
        """
        if isinstance(time, float):
            if time >= 0:
                if time < self.t_min:
                    return 0
                elif time <= END_TIME:
                    return np.interp(time, self.times, self.dtd)
                else:
                    return np.nan
            else:
                raise ValueError('Time must be positive.')
        else:
            raise TypeError('Parameter "time" must be a float. Got: %s' \
                            % type(time))

    @property
    def scheme(self):
        """
        Type: str
            Keyword denoting assumptions about the evolution of the
            double-degenerate system.
            Allowed values:

                - 'wide': there is a wide distribution of separation and total
                mass irrespective of the initial mass of the secondary
                - 'close': there is a tight correlation between the initial
                separation and the gravitational merge time
        """
        return self._scheme

    @scheme.setter
    def scheme(self, value):
        if isinstance(value, str):
            if value.lower() in ['wide', 'close']:
                self._scheme = value.lower()
            else:
                raise ValueError('DD scheme must be either "wide" or "close".')
        else:
            raise TypeError('Parameter "scheme" must be a string. Got: %s' % \
                            type(value))

    @property
    def efficiency(self):
        """
        Type: float
            The single-degenerate mass-transfer efficiency represented by
            epsilon in Equation 17 from Greggio 2005.
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
            raise TypeError('Parameter "efficiency" must be a Number. Got: %s'\
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
        value = check_if_instance(value, str, 'mlr')
        try:
            mlr_wrapper(1, model=value)
            self._mlr = value
        except:
            raise ValueError('Parameter "mlr" not in acceptable list.')

    @property
    def dt(self):
        """
        Type: float
            Integration timestep in Gyr
        """
        return self._dt

    @dt.setter
    def dt(self, value):
        if isinstance(value, Number):
            if value > 0:
                self._dt = float(value)
            else:
                raise ValueError('Integration timestep must be positive.')
        else:
            raise TypeError('Parameter "dt" must be a Number. Got: %s' \
                            % type(value))

    def integrate(self, time):
        """
        Integrate over all nuclear timescales up to the given time to obtain
        the distribution function of DD SNe Ia at that time.

        Parameters
        ----------
        time : float
            Time in Gyr since star formation event

        Returns
        -------
        float
            The value of the DTD at the given time
        """
        if time <= self.t_min or time >= self.t_max:
            return 0
        else:
            if self.scheme == 'wide':
                t_int_min = self.t_nuc_min
            else:
                t_int_min = self.asymptotic_nuclear_lifetime(time)
            t_int_max = min((self.t_nuc_max, time))
            integration_time = np.arange(t_int_min,
                                         t_int_max + self.dt,
                                         self.dt)
            # TODO this re-calculates the SD distribution every time, but it
            # doesn't change. figure out a way to speed this up
            integral = np.sum([
                self.single_degenerate_distribution(t_nuc) * \
                    self.f_merge(time, t_nuc, case=1) * \
                    self.dt \
                for t_nuc in integration_time])
            return integral

    def f_merge(self, time, t_nuc, case=1):
        """
        Distribution function of DD merger timescales for the 'wide' scheme.

        Parameters
        ----------
        time : float
            Time since starburst in Gyr
        t_nuc : float
            Nuclear timescale of the DD system in Gyr
        case : int [default: 1]
            Parameter controlling the form of the mass dependence for the
            'wide' scheme (see mass_dependence for more)

        Returns
        -------
        float
            Frequency of merger events
        """
        if self.scheme == 'wide':
            if t_nuc <= time - self.t_grav_min:
                m2 = mlr_wrapper(t_nuc, which='age', model=self.mlr)
                fW12 = self.mass_dependence(m2, case=case)
                return fW12 * (time - t_nuc) ** (-0.75 + 0.25 * self.beta)
            else:
                return 0
        else:
            if t_nuc <= time - self.t_grav_min:
                t_grav_max = self.maximum_gravitational_delay(t_nuc)
                t_diff = time - t_nuc
                return t_diff ** self.beta / (
                    t_grav_max ** (1 + self.beta) - \
                        self.t_grav_min ** (1 + self.beta))
            else:
                return 0

    def mass_dependence(self, m2, case=1):
        """
        Equivalent to the f_(1,2)^W term in Greggio 2005, defined in Equation 29.

        Parameters
        ----------
        m2 : float
            Mass of the secondary in solar masses
        case : int [default: 1]
            Parameter controlling the relationship between the secondary mass
            m2 and the total double-degenerate mass MDD. If case = 1, there
            is a tight correlation between m2 and MDD; if case = 2, a wide
            distribution of MDD is obtained from systems with the same m2
        """
        if case == 1:
            MDD = 1.4 + (m2 - 2)/6
            return MDD ** (0.75 + 0.75 * self.beta)
        elif case == 2:
            m_wd_min = minimum_wd_mass(m2, efficiency=self.efficiency)
            m2R = max((2, 2 + 10 * (m_wd_min - 0.6)))
            MDDn = max((1.4, m2R + 0.6))
            MDDx = m2R + 1.2
            exp = 1.75 + 0.75 * self.beta
            return MDDx ** exp - MDDn ** exp
        else:
            raise ValueError('Parameter "case" must be 1 or 2.')

    def asymptotic_nuclear_lifetime(self, time):
        """
        The lower integration limit for the 'close' DD scheme, represented in
        Greggio 2005 as t_(n,inf) and defined in equation 38.

        Parameters
        ----------
        time : float
            A delay time in Gyr

        Returns
        -------
        float
            The lower limit of integration in Gyr
        """
        if time < self.t_nuc_min + \
            self.maximum_gravitational_delay(self.t_nuc_min):
            return self.t_nuc_min
        else:
            # Solution to the equation t - t_nuc = t_(gw,x)(t_nuc)
            t_nuc_arr = np.arange(self.t_nuc_min, self.t_nuc_max, self.dt)
            lhs = time - t_nuc_arr
            rhs = np.array([self.maximum_gravitational_delay(t_nuc) \
                            for t_nuc in t_nuc_arr])
            soln = t_nuc_arr[np.argmin(np.abs(lhs - rhs))]
            return soln

    @staticmethod
    def maximum_gravitational_delay(t_nuc):
        """
        The upper envelope of the delay due to gravitational inspiral,
        represented in Greggio 2005 as tau_(gw,x) in Equation 30. The
        gravitational delay is proportional to the binary separation ^4.

        Parameters
        ----------
        t_nuc : float
            Nuclear lifetime of the secondary in Gyr.
        """
        log_t_nuc = m.log10(t_nuc) + 9
        log_t_grav_max = min((-16.66 + 3.17 * log_t_nuc,
                              6.02 + 0.52 * log_t_nuc))
        return 10 ** (log_t_grav_max - 9)


class greggio05_single:
    def __init__(self, m2_slope=-1.44, q_slope=1, efficiency=1.,
                 mlr='larson1974', imf='kroupa'):
        """
        Parameters
        ----------
        m2_slope : float [default: -1.44]
            Slope of the power-law derivative of secondary mass
        q_slope : float [default: 1]
            The power-law slope of the mass ratio distribution q = m2/m1
        efficiency : float [default: 0.5]
            Mass-transfer efficiency represented by epsilon in Equation 17
        mlr : str [default: 'larson1974']
            Which mass-lifetime relation to use
        imf : str [default: 'kroupa']
            Which IMF to use
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
        self._m2_slope = float(check_if_instance(value, Number, 'm2_slope'))

    @property
    def q_slope(self):
        """
        Type: float
            The power-law slope of the mass ratio distribution (q = m2/m1)
        """
        return self._q_slope

    @q_slope.setter
    def q_slope(self, value):
        self._q_slope = float(check_if_instance(value, Number, 'q_slope'))

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
        value = check_if_instance(value, Number, 'efficiency')
        if value >= 0 and value <= 1:
            self._efficiency = float(value)
        else:
            raise ValueError('Efficiency must be between 0 and 1.')

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
        value = check_if_instance(value, str, 'mlr')
        try:
            mlr_wrapper(1, model=value)
            self._mlr = value
        except:
            raise ValueError('Parameter "mlr" not in acceptable list.')

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
        value = check_if_instance(value, str, 'imf')
        if value.lower() in ['kroupa', 'salpeter']:
            self._imf = value
        else:
            raise ValueError('IMF must be either "kroupa" or "salpeter".')

    def secondary_mass_distribution(self, m2, m1_max=8, dm=0.01):
        """
        The distribution function of the secondaries in SN Ia progenitor
        systems (equation 16 in Greggio 2005). This assumes the mass of the
        primary follows a Kroupa or Salpeter distribution, and the mass ratio
        q = m2/m1 follows a power-law with a slope of 'q_slope'

        Parameters
        ----------
        m2 : float
            Mass of the secondary in solar masses
        m1_max : float [default: 8]
            Maximum primary mass in solar masses
        dm : float [default: 0.01]
            Mass integration step in solar masses

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

    def normalize(self, tmin=0.04, tmax=END_TIME, dt=1e-3):
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


def minimum_wd_mass(secondary_mass, efficiency=1.):
    """
    The minimum acceptable white dwarf mass for a successful SN Ia event.

    Parameters
    ----------
    secondary_mass : float
        Mass of the secondary in solar masses
    efficiency : float [default: 0.5]
        Efficiency of secondary mass transfer. In practice, any value
        above 0.3 will not affect the final DTD

    Returns
    -------
    float
        Minimum WD mass in solar masses
    """
    envelope_mass = secondary_mass - remnant_mass(secondary_mass)
    return 1.4 - efficiency * envelope_mass

def remnant_mass(initial_mass):
    """
    Relationship between a star's initial mass and remnant mass according
    to Case B Roche lobe overflow from Nelemans et al. 2001 (equation 18
    from Greggio 2005).

    Parameters
    ----------
    initial_mass : float
        Initial stellar mass in solar masses

    Returns
    -------
    float
        Remnant mass in solar masses
    """
    return max((0.3,
                0.3 + 0.1 * (initial_mass - 2),
                0.5 + 0.15 * (initial_mass - 4)))

def mlr_wrapper(qty, which='mass', model='larson1974', **kwargs):
    """
    Wrapper for the mass-lifetime relation which allows the selection of any
    MLR model available in VICE.

    Parameters
    ----------
    model : str
        Designator of a mass-lifetime relation in VICE
    All other arguments and keyword arguments are passed to vice.mlr.<model>

    Returns
    -------
    float
        Output of vice.mlr.<model>
    """
    return {
        'larson1974': vice.mlr.larson1974,
        'mm1989': vice.mlr.mm1989,
        'pm1993': vice.mlr.pm1993,
        'ka1997': vice.mlr.ka1997,
        'hpt2000': vice.mlr.hpt2000,
        'vincenzo2016': vice.mlr.vincenzo2016,
        'powerlaw': vice.mlr.powerlaw
    }[model.lower()](qty, which=which, **kwargs)

def check_if_instance(param, param_type, name):
    """
    Check if the given parameter is an instance of the given class.

    Parameters
    ----------
    param : any
        Parameter to check against the given class.
    param_type : type
        Type to check against parameter.
    name : str
        String representation of parameter name.

    Raises
    ------
    TypeError
        If 'param' is not an instance of 'param_type', 'param_type' is not an
        instance of 'type', or 'name' is not an instance of 'str'.

    Returns
    -------
    param : param_type
        Same as 'param' argument if it is of type 'param_type'.

    """
    if not isinstance(param_type, type):
        raise TypeError(
            'Parameter "%s" must be an instance of "%s". Got: "%s"' \
            % ('param_type', str(type), type(param_type)))
    if not isinstance(name, str):
        raise TypeError(
            'Parameter "%s" must be an instance of "%s". Got: "%s"' \
            % ('name', str(str), type(name)))
    if isinstance(param, param_type):
        return param
    else:
        raise TypeError(
            'Parameter "%s" must be an instance of "%s". Got: "%s"' \
            % (name, str(param_type), type(param)))
