"""
This script tests for bimodality in the [O/Fe] distribution of a VICE model
output by asking a simple question: is the [O/Fe] distribution function
within a given region of the Galaxy better fit by one Gaussian or two?
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from utils import import_apogee, apogee_region, get_bin_centers, kl_divergence
from utils import multioutput_to_pandas
from feh_distribution import apogee_mdf, vice_mdf
import paths

def main():
    # apogee_data = import_apogee(verbose=True)
    vice_stars = multioutput_to_pandas('diffusion/conroy22/exponential_timescale30')
    odf, bin_edges = vice_mdf(vice_stars, col='[o/fe]', galr_lim=(5, 7), 
                              absz_lim=(0.5, 1), xlim=(-0.15, 0.55), 
                              bin_width=0.01, smoothing=False)
    dx = bin_edges[1] - bin_edges[0]
    ofe_arr = get_bin_centers(bin_edges)
    # single Gaussian fit
    print('Fitting single Gaussian...')
    popt1, pcov1 = curve_fit(one_gaussian, ofe_arr, odf, p0=(0.1, 0.05))
    fit1 = one_gaussian(ofe_arr, *popt1)
    kld1 = kl_divergence(odf, fit1, dx)
    print('KL Divergence =', kld1)
    # double Gaussian fit
    print('Fitting double Gaussian...')
    popt2, pcov2 = curve_fit(two_gaussian, ofe_arr, odf, 
                             p0=(0.1, 0.05, 0.3, 0.05),
                             bounds=([-np.inf, -np.inf, 0.2, -np.inf],
                                     [0.2, 0.1, np.inf, 0.1]))
    fit2 = two_gaussian(ofe_arr, *popt2)
    kld2 = kl_divergence(odf, fit2, dx)
    print('KL Divergence =', kld2)
    # plot
    fig, ax = plt.subplots()
    ax.plot(ofe_arr, odf, label='VICE')
    ax.plot(ofe_arr, fit1, label='Single Gaussian fit')
    ax.plot(ofe_arr, fit2, label='Double Gaussian fit')
    ax.legend(loc='upper right')
    plt.savefig(paths.figures / 'bimodality_test.png', dpi=300)
    plt.close()
    
def one_gaussian(x, center, stdev):
    ng = NormalGaussian(center, stdev)
    return ng(x)

def two_gaussian(x, center1, stdev1, center2, stdev2):
    ng1 = NormalGaussian(center1, stdev1)
    ng2 = NormalGaussian(center2, stdev2)
    return 0.5 * ng1(x) + 0.5 * ng2(x)
    
class NormalGaussian:
    """
    A generic normalized Gaussian function of time.

    Attributes
    ----------
    center : float
        The location of the peak of the Gaussian function.
    stdev : float
        The standard deviation of the Gaussian function.
    norm : float
        The normalization coefficient of the Gaussian function.

    """
    def __init__(self, center=1, stdev=1):
        """
        Initialize the Gaussian.

        Parameters
        ----------
        center : float [default: 1]
            The location of the peak of the Gaussian function.
        stdev : float [default: 1]
            The standard deviation of the Gaussian function.
        """
        self.center = center
        self.stdev = stdev
        self.norm = 1 / (stdev * np.sqrt(2 * np.pi))

    def __call__(self, x):
        return self.norm * np.exp(-(x-self.center)**2 / (2*self.stdev**2))


if __name__ == '__main__':
    main()
