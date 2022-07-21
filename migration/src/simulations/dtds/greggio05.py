"""
This file contains a prescription for the DTD of SNe Ia following the double-
degenerate solution of Greggio 2005.
"""

import math as m
import vice

def main():
    pass

class greggio05:
    def __init__(self, scheme='wide', t_nuc_min=0.04, t_nuc_max=1,
                 t_grav_min=0.001, dt=0.01):
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
        """
        self.t_nuc_min = t_nuc_min
        minimum_delay_total = t_nuc_min + t_grav_min

        if scheme.lower() == 'wide':
            t_min = t_nuc_min
        elif scheme == 'close':
            # t_min = asymptotic_nuclear_lifetime(t)
            pass
        else:
            raise ValueError('Parameter "scheme" must be "wide" or "close".')

    def __call__(self, time):
        pass

    def asymptotic_nuclear_lifetime(self, t):
        """
        The lower integration limit for the 'close' DD scheme, represented in
        Greggio 2005 as t_(n,inf).

        Parameters
        ----------
        t : float
            A delay time in Gyr

        Returns
        -------
        float
            The lower limit of integration in Gyr
        """
        if t < self.t_nuc_min + maximum_grativational_delay(t_nuc_min):
            return self.t_nuc_min
        else:
            pass

    def maximum_gravitational_delay(t_nuc):
        """
        The upper envelope of the delay due to gravitational inspiral,
        represented in Greggio 2005 as tau_(gw,x). The gravitational delay
        is proportional to the binary separation ^4.

        Parameters
        ----------
        t_nuc : float
            Nuclear lifetime of the secondary in Gyr.
        """
        log_t_nuc = m.log10(t_nuc) + 9
        log_t_grav_max = min((-16.66 + 3.17 * log_t_nuc, 6.02 + 0.52 * log_t_nuc))
        return 10 ** (log_t_grav_max - 9)

    def secondary_mass_distribution(m2, m2_max=8):
        pass

if __name__ == '__main__':
    main()
