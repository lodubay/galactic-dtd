"""
This script plots the distribution of stellar ages as predicted by VICE.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vice
import paths

def main():
    pass

def f_survive(age, mlr='larson1974', imf='kroupa', m_lower=0.08, m_upper=100,
              dm=0.01):
    """
    Calculate the surviving mass fraction of a stellar population with a given
    age.

    Parameters
    ----------
    age : float
        Age of the stellar population in Gyr
    mlr : str
        Mass-lifetime relation (MLR) to use. The default is 'larson1974'.
    imf : str
        Which IMF to use. Options are 'kroupa' or 'salpeter'.
    m_lower : float
        Lower limit of stellar mass. The default is 0.08 solar masses.
    m_upper : float
        Upper limit of stellar mass. The default is 100 solar masses.
    dm : float
        Integration step size in solar masses. The default is 0.01.

    Returns
    -------
    float
        Stellar population surviving mass fraction.
    """
    mlr_select = {
        'larson1974': vice.mlr.larson1974,
        'mm1989': vice.mm1989,
        'pm1993': vice.mlr.pm1993,
        'ka1997': vice.mlr.ka1997,
        'hpt2000': vice.mlr.hpt2000,
        'vincenzo2016': vice.mlr.vincenzo2016,
        'powerlaw': vice.mlr.powerlaw
    }
    if mlr in mlr_select.keys():
        # Mass of a star with a lifetime equal to age
        mass = mlr_select[mlr](age, which='age')
        m_arr = np.array(m_lower, min((m_upper, mass)), dm)
        normal_imf = NormalIMF(which=imf, m_lower=m_lower, m_upper=m_upper,
                               dm=dm)
        f_survive = sum([normal_imf(m) for m in m_arr])
        return f_survive
    else:
        raise ValueError('MLR must be in acceptable list.')


class NormalIMF:
    """
    A normalized initial mass function (IMF).
    """
    def __init__(self, which='kroupa', **kwargs):
        """
        Initialize the IMF.

        Parameters
        ----------
        which : string, optional
            Which version of the IMF to use. Must be one of 'salpeter' or
            'kroupa'. The default is 'kroupa'
        Other keyword arguments are passed to self.normalize
        """
        select = {
            'salpeter': vice.imf.salpeter,
            'kroupa': vice.imf.kroupa
        }
        if which in select.keys():
            self._imf = select[which]
            self.norm = self.normalize(**kwargs)
        else:
            raise ValueError('IMF must be either "kroupa" or "salpeter".')

    def __call__(self, mass):
        """
        Calculate the normalized IMF at a given stellar mass.

        Parameters
        ----------
        mass : float
            Stellar mass in solar masses.

        Returns
        -------
        float
            The normalized value of the IMF at that stellar mass.
        """
        return self.norm * self._imf(mass)

    def normalize(self, m_lower=0.08, m_upper=100, dm=0.01):
        """
        Calculate the normalization coefficient of the IMF.

        Parameters
        ----------
        m_lower : float, optional
            Lower limit of integration in solar masses. The default is 0.08.
        m_upper : TYPE, optional
            Upper limit of integration in solar masses. The default is 100.
        dm : TYPE, optional
            Integration step in solar masses. The default is 0.01.

        Returns
        -------
        float
            Normalization coefficient of the IMF.
        """
        m_arr = np.arange(m_lower, m_upper, dm)
        integral = sum([self._imf(m) * dm for m in m_arr])
        return 1 / integral


if __name__ == '__main__':
    main()
