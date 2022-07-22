"""
This file contains a prescription for the DTD of SNe Ia following the double-
degenerate solution of Greggio 2005.
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import vice

def main():
    tarr = np.logspace(np.log10(0.03), np.log10(13.2), num=50, base=10)
    dtd_wide = greggio05(scheme='wide')
    RIa_wide = [dtd_wide(t) for t in tarr]
    plt.plot(tarr * 1e9, RIa_wide, label='wide kroupa')
    # dtd_wide_salpeter = greggio05(scheme='wide', imf='salpeter')
    # RIa_wide_salpeter = [dtd_wide_salpeter(t) for t in tarr]
    # plt.plot(tarr * 1e9, RIa_wide_salpeter, ls='--', label='wide salpeter')
    dtd_close = greggio05(scheme='close')
    RIa_close = [dtd_close(t) for t in tarr]
    plt.plot(tarr * 1e9, RIa_close, label='close kroupa')
    # dtd_close_salpeter = greggio05(scheme='close', imf='salpeter')
    # RIa_close_salpeter = [dtd_close_salpeter(t) for t in tarr]
    # plt.plot(tarr * 1e9, RIa_close_salpeter, ls='--', label='close salpeter')
    plt.xlabel('Time [yr]')
    plt.xscale('log')
    plt.yscale('log')
    # plt.ylim((3e-3, 3))
    plt.legend()
    plt.show()

class greggio05:
    def __init__(self, scheme='wide', t_nuc_min=0.04, t_nuc_max=1,
                 t_grav_min=0.001, beta_sep=0, beta_grav=-0.75, dt=0.001,
                 mlr='larson1974', imf='kroupa'):
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
        mlr : str [default: 'larson1974']
            Which model of mass-lifetime relation to use. Must be one available
            in VICE
        imf : str [default: 'kroupa']
            Which model of IMF to use. Must be one available in VICE
        """
        self.t_nuc_min = t_nuc_min
        self.t_nuc_max = t_nuc_max
        self.t_grav_min = t_grav_min
        self.dt = dt
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
        Parameters
        ----------
        time : float
            Time after starburst in Gyr

        Returns
        -------
        float
            Value of the distribution function of DD SNe Ia
        """
        print(time)
        if isinstance(time, float):
            if time <= self.t_min or time >= self.t_max:
                return 0
            else:
                if self.scheme == 'wide':
                    t_int_min = self.t_nuc_min
                    f_merge = self.f_merge_wide
                else:
                    t_int_min = self.asymptotic_nuclear_lifetime(time)
                    f_merge = self.f_merge_close
                t_int_max = min((self.t_nuc_max, time))
                integration_time = np.arange(t_int_min,
                                             t_int_max + self.dt,
                                             self.dt)
                integral = sum([
                    single_degenerate_distribution(
                            t_nuc, mlr=self.mlr, imf=self.imf, slope=-0.7) * \
                        f_merge(time, t_nuc, beta=self.beta) * self.dt \
                        for t_nuc in integration_time])
                return integral
        else:
            raise TypeError('Parameter "time" must be a float. Got: %s' \
                            % type(time))

    def f_merge_wide(self, time, t_nuc, case=1, beta=0, mlr='larson1974'):
        """
        Distribution function of DD merger timescales for the 'wide' scheme.

        Parameters
        ----------
        time : float
            Time since starburst in Gyr
        t_nuc : float
            Nuclear timescale of the DD system in Gyr
        case : int [default: 1]
            Parameter controlling the form of the mass dependence
            (see mass_dependence for more)
        beta : float [default: 0]
            Slope of the power-law distribution of separation
        mlr : str [default: 'larson1974']
            Which model of mass-lifetime relation to use

        Returns
        -------
        float
            Frequency of merger events
        """
        if t_nuc <= time - self.t_grav_min:
            m2 = mlr_wrapper(t_nuc, which='age', model=mlr)
            fW12 = self.mass_dependence(m2, case=case)
            return fW12 * (time - t_nuc) ** (-0.75 + 0.25 * beta)
        else:
            return 0

    def f_merge_close(self, time, t_nuc, beta=-0.75):
        """
        Distribution function of DD merger timescales for the 'close' scheme.

        Parameters
        ----------
        time : float
            Time since starburst in Gyr
        t_nuc : float
            Nuclear timescale of the DD system in Gyr
        beta : float [default: -0.75]
            Slope of the power-law distribution of gravitational delay

        Returns
        -------
        float
            Frequency of merger events
        """
        if t_nuc <= time - self.t_grav_min:
            t_grav_max = self.maximum_gravitational_delay(t_nuc)
            t_diff = time - t_nuc
            return t_diff ** beta / (
                t_grav_max ** (1 + beta) - self.t_grav_min ** (1 + beta))
        else:
            return 0

    def mass_dependence(self, m2, case=1, beta=0):#, minimum_wd_mass=0.7):
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
        beta : float [default: 0]
            Slope of the power-law distribution of separation
        """
        if case == 1:
            MDD = 1.4 + (m2 - 2)/6
            return MDD ** (0.75 + 0.75 * beta)
        elif case == 2:
            m2R = max((2, 2 + 10 * (minimum_wd_mass(m2) - 0.6)))
            MDDn = max((1.4, m2R + 0.6))
            MDDx = m2R + 1.2
            exp = 1.75 + 0.75 * beta
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

def single_degenerate_distribution(time, mlr='larson1974', imf='kroupa',
                                   slope=-1.44):
    """
    The DTD of the Greggio 2005 single-degenerate case.

    Parameters
    ----------
    time : float
        Time after starburst in Gyr
    mlr : str [default: 'larson1974']
        Which mass-lifetime relation to use
    slope : float [default: -1.44]
        Slope of the power-law derivative of secondary mass

    Returns
    -------
    float
        Frequency of SD events
    """
    secondary_mass = mlr_wrapper(time, which='age', model=mlr)
    # log_m2 = 0.0471 * (np.log10(time * 1e9))**2 - 1.2 * np.log10(time * 1e9) + 7.3
    # secondary_mass = 10 ** log_m2
    return secondary_mass_distribution(secondary_mass, imf=imf) * time ** slope

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

def secondary_mass_distribution(m2, m1_max=8, imf='kroupa',
                                q_slope=1,dm=0.01):
    # minimum_wd_mass=0.7, dm=0.01):
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
    m1_min = max((m2, 2, 2 + 10 * (minimum_wd_mass(m2) - 0.6)))
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

if __name__ == '__main__':
    main()
