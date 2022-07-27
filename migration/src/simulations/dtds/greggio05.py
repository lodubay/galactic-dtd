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
                 t_grav_min=0.001, beta_sep=0., beta_grav=-0.75, dt=0.001,
                 mlr='larson1974', imf='kroupa', sd_slope=-1.44, nsamples=100,
                 q_slope=1, efficiency=0.5):
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
        t_grav_min : float [default: 0.001]
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
        dt : float [default: 0.001]
            Integration timestep in Gyr
        mlr : str [default: 'larson1974']
            Which model of mass-lifetime relation to use. Must be one available
            in VICE
        imf : str [default: 'kroupa']
            Which model of IMF to use. Must be one available in VICE
        sd_slope : float [default: -1.44]
            Power-law slope of the secondary mass evolution function, which
            goes into the single-degenerate distribution function
        nsamples : int [default: 100]
            Number of log-spaced points at which to pre-compute the DTD
        q_slope : float [default: 1]
            The power-law slope of the mass ratio distribution q = m2/m1
        efficiency : float [default: 0.5]
            Mass-transfer efficiency represented by epsilon in Equation 17
        """
        self.t_nuc_min = check_if_instance(t_nuc_min, float, 't_nuc_min')
        self.t_nuc_max = check_if_instance(t_nuc_max, float, 't_nuc_max')
        self.t_grav_min = check_if_instance(t_grav_min, float, 't_grav_min')
        if dt > 0:
            self.dt = check_if_instance(dt, float, 'dt')
        else:
            raise ValueError('Integration timestep must be positive.')
        # Total minimum and maximum delay times
        self.t_min = t_nuc_min + t_grav_min
        self.t_max = t_nuc_max + self.maximum_gravitational_delay(t_nuc_max)

        if scheme.lower() == 'wide':
            self.beta = beta_sep
        elif scheme == 'close':
            # t_min = asymptotic_nuclear_lifetime(t)
            self.beta = beta_grav
        else:
            raise ValueError('Parameter "scheme" must be "wide" or "close".')
        self.scheme = scheme.lower()

        if efficiency >= 0 and efficiency <= 1:
            self.efficiency = check_if_instance(efficiency, Number, 'efficiency')
        else:
            raise ValueError('Efficiency must be a number between 0 and 1.')

        if isinstance(mlr, str):
            try:
                mlr_wrapper(1, model=mlr)
                self.mlr = mlr
            except:
                raise ValueError('Parameter "mlr" not in acceptable list.')
        else:
            raise TypeError('Parameter "mlr" must be a string. Got: %s' \
                            % type(mlr))

        self.single_degenerate_distribution = greggio05_single(
            m2_slope=sd_slope, q_slope=q_slope, efficiency=efficiency,
            mlr=mlr, imf=imf
        )

        # Pre-calculate the DTD to save on compute time
        print('Computing Greggio 2005 DTD...')
        self.times = np.logspace(np.log10(self.t_min),
                                 np.log10(min((END_TIME, self.t_max))),
                                 num=nsamples, endpoint=True)
        self.dtd = np.zeros(self.times.shape)
        for i, t in enumerate(tqdm(self.times)):
            self.dtd[i] = self.pre_calculate(t)


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
            Interpolated value of the distribution function of DD SNe Ia
        """
        if isinstance(time, float):
            return np.interp(time, self.times, self.dtd)
        else:
            raise TypeError('Parameter "time" must be a float. Got: %s' \
                            % type(time))

    def pre_calculate(self, time, t_min=0):
        """
        Integrate over all nuclear timescales up to the given time to obtain
        the distribution function of DD SNe Ia at that time.

        Parameters
        ----------
        time : float
            Time in Gyr since star formation event
        t_min : float [default: 0]
            Minimum limit of integration, which will override the default if
            it is greater than the automatically calculated limit

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
            integral = sum([
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
    def __init__(self, m2_slope=-1.44, q_slope=1, efficiency=0.5,
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
        if isinstance(m2_slope, Number):
            self.m2_slope = m2_slope
        else:
            raise TypeError('Parameter "m2_slope" must be a number. Got: %s' \
                            % type(m2_slope))
        if isinstance(q_slope, Number):
            self.q_slope = q_slope
        else:
            raise TypeError('Parameter "q_slope" must be a number. Got: %s' \
                            % type(q_slope))
        if isinstance(efficiency, Number):
            if efficiency >= 0 and efficiency <= 1:
                self.efficiency = efficiency
            else:
                raise ValueError('Efficiency must be between 0 and 1.')
        else:
            raise TypeError('Parameter "efficiency" must be a number. Got: %s'\
                            % type(efficiency))
        if isinstance(mlr, str):
            try:
                mlr_wrapper(1, model=mlr)
                self.mlr = mlr
            except:
                raise ValueError('Parameter "mlr" not in acceptable list.')
        else:
            raise TypeError('Parameter "mlr" must be a string. Got: %s' \
                            % type(mlr))
        if isinstance(imf, str):
            if imf.lower() in ['kroupa', 'salpeter']:
                self.imf = imf
            else:
                raise ValueError('Parameter "imf" must be one of "kroupa" or ' + \
                                 '"salpeter".')
        else:
            raise TypeError('Parameter "imf" must be a string. Got: %s' \
                            % type(imf))

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
        f_m2 = self.secondary_mass_distribution(secondary_mass, imf=self.imf,
                                                q_slope=self.q_slope,
                                                efficiency=self.efficiency)
        return f_m2 * time ** self.m2_slope

    @staticmethod
    def secondary_mass_distribution(m2, m1_max=8, imf='kroupa', q_slope=1,
                                    dm=0.01, efficiency=0.5):
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
        imf : str [default: 'kroupa']
            Which IMF to assume, either 'kroupa' or 'salpeter'
        q_slope : float [default: 1]
            The power-law slope of the mass ratio distribution q = m2/m1
        dm : float [default: 0.01]
            Mass integration step in solar masses

        Returns
        -------
        float
            The relative frequency of the given secondary mass
        """
        # Lower limit of integration over primary mass
        m_wd_min = minimum_wd_mass(m2, efficiency=efficiency)
        m1_min = max((m2, 2, 2 + 10 * (m_wd_min - 0.6)))
        primary_masses = np.arange(m1_min, m1_max + dm, dm)
        imf_func = {
            'kroupa': vice.imf.kroupa,
            'salpeter': vice.imf.salpeter
        }[imf]
        integral = sum([m2**q_slope * m1**q_slope * imf_func(m1) * dm \
                        for m1 in primary_masses])
        return integral


def minimum_wd_mass(secondary_mass, efficiency=0.5):
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
    }[model](qty, which=which, **kwargs)

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
