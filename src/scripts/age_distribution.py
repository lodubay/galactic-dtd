"""
This script plots the distribution of stellar ages as predicted by VICE.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vice
import paths
from utils import multioutput_to_pandas
from _globals import DT

def main():
    output = 'diffusion/insideout/powerlaw'
    print('import stars')
    stars = multioutput_to_pandas(output)
    print('age distribution')
    ages, n_stars = age_distribution(stars)
    plt.plot(ages, n_stars)
    plt.xlabel('Age [Gyr]')
    plt.ylabel('dN/dAge')
    plt.show()
    # age = 0.01
    # imf = NormalIMF(m_upper=vice.mlr.larson1974(age, which='age'))
    # mass = 5000
    # print(imf.weighted_mean())


def age_distribution(stars, range=None, bin_width=1, **kwargs):
    if (bin_width / DT) % 1 == 0:
        if not range:
            range = (stars['age'].min(), stars['age'].max())
        bins = np.arange(range[0], range[1] + bin_width, bin_width)
        # Sum stellar mass in each timestep
        mass_total = stars.groupby(['age']).sum()['mass']
        ages = np.array(mass_total.index)
        mass_total = np.array(mass_total)
        # Calculate remaining stellar mass today
        mass_remaining = mass_total * (1 - np.array(
            [vice.cumulative_return_fraction(age) for age in ages]))
        # Average mass of a star of that particular age
        mass_average = np.array([mean_stellar_mass(age) for age in ages])
        # Number of stars in each age bin
        n_stars = np.around(mass_remaining / mass_average)
        n_stars /= DT
        return ages, n_stars
    else:
        raise ValueError('Bin width must be a multiple of the simulation timestep.')


def mean_stellar_mass(age, imf=vice.imf.kroupa, mlr=vice.mlr.larson1974,
                      m_lower=0.08, m_upper=100, dm=0.01):
    """
    Calculate the mean mass of a stellar population of a given age.

    Parameters
    ----------
    age : float
        Stellar age in Gyr
    imf : <function>, optional
        Initial mass function which takes mass in solar masses as an argument.
        The default is vice.imf.kroupa
    mlr : <function>, optional
        Mass-lifetime relation which takes age in Gyr as an argument. The
        default is vice.mlr.larson1974
    m_lower : float, optional
        Lower mass limit on IMF in solar masses. The default is 0.08
    m_upper : float, optional
        Upper mass limit on IMF in solar masses. The default is 100
    dm : float, optional
        IMF integration step in solar masses. The default is 0.01

    Returns
    -------
    float
        Mean mass of stars with lifetime greater than or equal to the given age
        weighted by the IMF
    """
    m_max = min((mlr(age, which='age'), m_upper))
    masses = np.arange(m_lower, m_max + dm, dm)
    dndm = np.array([imf(m) for m in masses])
    weighted_mean = np.average(masses, weights=dndm)
    return weighted_mean


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
    def __init__(self, which='kroupa', m_lower=0.08, m_upper=100, dm=0.01):
        """
        Initialize the IMF.

        Parameters
        ----------
        which : string, optional
            Which version of the IMF to use. Must be one of 'salpeter' or
            'kroupa'. The default is 'kroupa'
            m_lower : float, optional
                Lower limit of integration in solar masses. The default is 0.08.
            m_upper : TYPE, optional
                Upper limit of integration in solar masses. The default is 100.
            dm : TYPE, optional
                Integration step in solar masses. The default is 0.01.
        """
        select = {
            'salpeter': vice.imf.salpeter,
            'kroupa': vice.imf.kroupa
        }
        if which in select.keys():
            self.dm = dm
            self.masses = np.arange(m_lower, m_upper + dm, dm)
            self._imf = select[which]
            self.norm = 1
            self.norm = 1 / self.integrate()
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

    def integrate(self):
        """
        float
            The integral of the IMF
        """
        integral = sum([self.__call__(m) * self.dm for m in self.masses])
        return integral

    def weighted_mean(self):
        """
        Calculate the average stellar mass of the IMF.
        """
        weights = np.array([self.__call__(m) for m in self.masses])
        weighted_mean = np.average(self.masses, weights=weights)
        return weighted_mean


if __name__ == '__main__':
    main()
