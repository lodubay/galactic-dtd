"""
This file contains a prescription for the DTD of SNe Ia following the double-
degenerate solution of Greggio 2005.
"""

from numbers import Number
import math as m
import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit
import vice
import matplotlib.pyplot as plt
from ..._globals import END_TIME

# The minimum nuclear lifetime of the secondary in Gyr (8 solar masses)
T_NUC_MIN = 0.04
# The maximum nuclear lifetime of the secondary in Gyr (2 solar masses)
T_NUC_MAX = 1.
# The minimum gravitational inspiral delay in Gyr
T_GRAV_MIN = 1e-3

def main():
    fit = greggio05_analytic.fit_to_model('wide', dt=1e-3, nsamples=100)
    print('Slope 1 = %s\nSlope 2 = %s' % (fit.slope1, fit.slope2))
    test = greggio05_analytic()
    tarr = np.arange(0.04, 13.2, 0.001)
    plt.plot(tarr * 1e9, [test(t) for t in tarr])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

class greggio05_analytic:
    """
    An analytic approximation of the double-degenerate DTD from Greggio 2005
    """
    def __init__(self, tsplit=1., slope1=-0.2, slope2=-0.9,
                 rise_strength=1.7, rise_timescale=0.09,
                 tail_strength=0.8, tail_timescale=0.1, normalize=True):
        self.tsplit = tsplit
        self.slope1 = slope1
        self.slope2 = slope2
        self.rise = self.exponential(rise_strength, rise_timescale)
        self.tail = self.exponential(-tail_strength, tail_timescale,
                                     offset=self.tsplit)
        self.norm = 1
        if normalize:
            self.norm = self.normalize()

    def __call__(self, time):
        if isinstance(time, Number):
            if time < 0:
                raise ValueError('Time must be non-negative.')
            elif time <= self.tsplit:
                return self.norm * self.part1(time)
            else:
                return self.norm * self.part2(time) * \
                    (self.part1(self.tsplit) / self.part2(self.tsplit))
        else:
            raise TypeError('Parameter "time" must be numeric. Got: %s' \
                            % type(time))

    @classmethod
    def fit_to_model(cls, scheme, tmin=0.04, tmax=END_TIME, tstep=0.01,
                     **kwargs):
        """
        Fit the analytic function to the greggio05_double model from scratch.

        Warning! Takes a long time.
        """
        model = greggio05_double(scheme, **kwargs)
        tarr = np.arange(tmin, tmax, tstep)
        yarr = np.array([model(t) for t in tarr])
        return cls.fit_to_data(tarr, yarr)

    @classmethod
    def fit_to_data(cls, tarr, yarr):
        """
        Fit the analytic function to data from a pre-generated DTD.
        """
        print('Fitting analytic model...')
        popt, pcov = curve_fit(analytic_wrapper, tarr, yarr,
                               p0=(0, -1, 1, 0.1, 1, 0.1, 1e-9))
        return cls(slope1=popt[0], slope2=popt[1],
                   rise_strength=popt[2], rise_timescale=popt[3],
                   tail_strength=popt[4], tail_timescale=popt[5])

    def normalize(self, tmin=0.04, tmax=END_TIME, dt=0.001):
        tarr = np.arange(tmin, tmax+dt, dt)
        integral = np.sum([self.__call__(t) * dt * 1e9 for t in tarr])
        return 1 / integral

    def part1(self, time):
        return self.rise(time) * time ** self.slope1

    def part2(self, time):
        return self.tail(time) * time ** self.slope2

    @staticmethod
    def exponential(strength, timescale, offset=0):
        """
        Generate an exponential scaling function of time.

        Parameters
        ----------
        strength : float
            The strength of the exponential scaling. If positive, the function
            increases over time; if negative, it decreases over time
        timecale : float
            The exponential timescale in the same units of time
        offset : float, optional
            Horizontal (temporal) offset of the exponential. The default is 0

        Returns
        -------
        <function>
            A function of time
        """
        return lambda t: (1 - strength * m.exp(-(t - offset) / timescale))

def analytic_wrapper(time, slope1, slope2, rise_strength, rise_timescale,
                     tail_strength, tail_timescale, norm):
    cl = greggio05_analytic(slope1=slope1, slope2=slope2,
                            rise_strength=rise_strength,
                            rise_timescale=rise_timescale,
                            tail_strength=tail_strength,
                            tail_timescale=tail_timescale,
                            normalize=False)
    if isinstance(time, np.ndarray):
        return np.array([norm * cl(t) for t in time])
    else:
        return norm * cl(time)


class greggio05_double:
    """
    An numerical implementation of the double-degenerate DTD from Greggio 2005

    Attributes
    ----------
    scheme : str
        Keyword denoting assumptions about the evolution of a double-degenerate
        system, either 'wide' or 'close'
    beta : float
        The exponent of the power-law distribution of merger timescales
    efficiency : float
        The efficiency of secondary mass transfer
    dt : float
        The integration timestep in Gyr
    mlr : str
        Which of the VICE mass-lifetime relations to assume
    t_min : float
        Total minimum delay time in Gyr
    t_max : float
        Total maximum delay time in Gyr
    single_degenerate_distribution : greggio05_single
        Instance of the single-degenerate DTD from Greggio 2005
    times : numpy.ndarray
        Array of logarithmically-spaced times at which the DTD is calculated
    dtd : numpy.ndarray
        Array of pre-calculated normalized values of the DTD

    Methods
    -------
    integrate(time)
        Calculate the value of the DTD at the given time
    f_merge(self, time, t_nuc, case=1)
        The distribution function of DD merger timescales
    mass_dependence(self, m2, case=1)
        The term describing the dependence on the mass of the DD system
    asymptotic_nuclear_lifetime(self, time)
        The lower integration limit for the 'close' DD scheme
    maximum_gravitational_delay(t_nuc)
        Static method to calculate the maximum gravitational delay time
    """

    def __init__(self, scheme, beta_sep=0., beta_grav=-0.75, efficiency=1.,
                 dt=1e-3, nsamples=100, mlr='larson1974', **kwargs):
        """
        Parameters
        ----------
        scheme : str
            Keyword denoting assumptions about the evolution of a DD system,
            either 'wide' or 'close'.
        beta_sep : float, optional
            The exponent of the power law distribution adopted for the final
            separations in the 'wide' scheme; a flat distribution of
            of separations corresponds to beta_sep = 0. The default is 0.
        beta_grav : float, optional
            The exponent of the power law distribution adopted for the
            gravitational delays in the 'close' scheme; a flat distribution
            of gravitational delays corresponds to beta_grav = -0.75. The
            default is -0.75.
        efficiency : float, optional
            The mass-transfer efficiency represented by epsilon in Equation 17.
            The default is 1.
        dt : float, optional
            Integration timestep in Gyr. The default is 1e-3.
        nsamples : int, optional
            Number of logarithmically-spaced points at which to pre-compute the
            DTD. The default is 100.
        mlr : str, optional
            Which model of mass-lifetime relation to use; must be one available
            in VICE. The default is 'larson1974'.
        kwargs : dict
            Keyword arguments passed to greggio05_single
        """

        self.scheme = scheme
        if scheme == 'wide':
            self.beta = beta_sep
        else:
            self.beta = beta_grav
        self.efficiency = efficiency
        self.dt = dt
        self.mlr = mlr

        # Total minimum and maximum delay times
        self.t_min = T_NUC_MIN + T_GRAV_MIN
        self.t_max = T_NUC_MAX + self.maximum_gravitational_delay(T_NUC_MAX)

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

        Raises
        ------
        ValueError
            If the time passed as a parameter is not positive.
        TypeError
            If the time passed as a parameter is not a float.
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
    def beta(self):
        """
        Type: float
            The exponent of the power-law distribution of merger timescales.
            If scheme == 'wide', this is assumed to be the distribution of
            final separations; if scheme == 'close', this is assumed to be the
            distribution of gravitational delay times.
        """
        return self._beta

    @beta.setter
    def beta(self, value):
        if isinstance(value, Number):
            self._beta = float(value)
        else:
            raise TypeError('Parameter "beta" must be numeric. Got: %s' \
                            % type(value))

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
            raise TypeError('Parameter "efficiency" must be numeric. Got: %s'\
                            % type(value))

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
            raise TypeError('Parameter "dt" must be numeric. Got: %s' \
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

    def integrate(self, time):
        """
        Calculate the value of the DTD at the given time

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
                t_int_min = T_NUC_MIN
            else:
                t_int_min = self.asymptotic_nuclear_lifetime(time)
            t_int_max = min((T_NUC_MAX, time))
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
        The distribution function of DD merger timescales.

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
            if t_nuc <= time - T_GRAV_MIN:
                m2 = mlr_wrapper(t_nuc, which='age', model=self.mlr)
                fW12 = self.mass_dependence(m2, case=case)
                return fW12 * (time - t_nuc) ** (-0.75 + 0.25 * self.beta)
            else:
                return 0
        else:
            if t_nuc <= time - T_GRAV_MIN:
                t_grav_max = self.maximum_gravitational_delay(t_nuc)
                t_diff = time - t_nuc
                return t_diff ** self.beta / (
                    t_grav_max ** (1 + self.beta) - \
                        T_GRAV_MIN ** (1 + self.beta))
            else:
                return 0

    def mass_dependence(self, m2, case=1):
        """
        The term describing the dependence on the mass of the DD system

        Equivalent to the f_(1,2)^W term in Greggio 2005, defined in
        Equation 29.

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
        The lower integration limit for the 'close' DD scheme

        Represented in Greggio 2005 as t_(n,inf) and defined in equation 38.

        Parameters
        ----------
        time : float
            A delay time in Gyr

        Returns
        -------
        float
            The lower limit of integration in Gyr
        """
        if time < T_NUC_MIN + \
            self.maximum_gravitational_delay(T_NUC_MIN):
            return T_NUC_MIN
        else:
            # Solution to the equation t - t_nuc = t_(gw,x)(t_nuc)
            t_nuc_arr = np.arange(T_NUC_MIN, T_NUC_MAX, self.dt)
            lhs = time - t_nuc_arr
            rhs = np.array([self.maximum_gravitational_delay(t_nuc) \
                            for t_nuc in t_nuc_arr])
            soln = t_nuc_arr[np.argmin(np.abs(lhs - rhs))]
            return soln

    @staticmethod
    def maximum_gravitational_delay(t_nuc):
        """
        Calculate the maximum gravitational delay time.

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
    efficiency : float, optional
        Efficiency of secondary mass transfer. In practice, any value
        above 0.3 will not affect the final DTD. The default is 1.

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

if __name__ == '__main__':
    main()
