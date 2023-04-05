"""
This script tests for bimodality in the [O/Fe] distribution of a VICE model
output by fitting a Gaussian to the distribution of an inner and outer region
of the Galaxy and calculating the difference between them.
"""

from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm, skewnorm, chisquare
import matplotlib.pyplot as plt
from utils import multioutput_to_pandas, get_bin_centers, kl_divergence
from feh_distribution import vice_mdf
import paths

OFE_LIM = (-0.15, 0.55)
BIN_WIDTH = 0.01
SIG = 3

SFHs = ['insideout', 
        'lateburst', 
        'conroy22', 
        'conroy22_JW20yields', 
        'twoinfall']
DTDs = ['powerlaw_slope11', 
        'powerlaw_slope14', 
        'exponential_timescale15',
        'exponential_timescale30',
        'plateau_width300_slope11',
        'plateau_width1000_slope11',
        'prompt_peak050_stdev015_timescale30',
        'triple_delay040']

def main():
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['DeltaChi2', 'KLDiv', 'DeltaMu', 'LowAlphaCoeff'],
    )
    
    with tqdm(total=len(SFHs) * len(DTDs)) as t:
        for evolution in SFHs:
            for RIa in DTDs:
                output_name = '/'.join(('diffusion', evolution, RIa))
                summary_table.loc[RIa, evolution] = bimodality(output_name)
                t.update()
    
    print(summary_table)
    summary_table.to_csv(paths.figures / 'bimodality/skewnorm/bimodality.csv')
    

def bimodality(output_name, diagnostic=True, galr_lim=(5, 7), absz_lim=(0.5, 1)):
    """
    Determine if the [O/Fe] distribution of the given VICE output is bimodal.
    
    This function fits a single skewed Gaussian and a model with two skewed
    Gaussians to the [O/Fe] distribution in a galactic region known to show
    alpha-element bimodality. If the two-Gaussian fit has a lower chi-square
    statistic than the one-Gaussian fit, the distribution is bimodal.
    
    Parameters
    ----------
    output_name : str
        Name of VICE multizone output, e.g. 
        'diffusion/insideout/powerlaw_slope11'
    diagnostic : bool, optional
        If True, make a diagnostic plot illustrating the inner and outer ODFs
        and their respective fits. The default is True.
    galr_lim : tuple, optional
        Bounds on the galactocentric radius in kpc. The default is (5, 7).
    absz_lim : tuple, optional
        Bounds on the distance from the midplane in kpc. The default is (0.5, 1).
    
    Returns
    -------
    delta_chisq : float
        Difference in the chi-square statistic between the one- and two-Gaussian
        fits. If this is positive, the distribution is bimodal.
    """
    vice_stars = multioutput_to_pandas(output_name)
    # galactic region ODF
    odf, bin_edges = vice_mdf(vice_stars, col='[o/fe]', 
                              galr_lim=galr_lim, absz_lim=absz_lim, 
                              xlim=OFE_LIM, bin_width=BIN_WIDTH, 
                              smoothing=False)
    dx = bin_edges[1] - bin_edges[0]
    ofe_arr = get_bin_centers(bin_edges)
    
    # one-Gaussian fit
    try:
        popt1, pcov1 = curve_fit(skewnorm.pdf, ofe_arr, odf, 
                                 p0=(0, 0.1, 0.05),
                                 bounds=([-np.inf, -np.inf, 0.], 
                                         [np.inf, np.inf, np.inf]))
        fit1 = skewnorm.pdf(ofe_arr, *popt1)
    # if skew-normal fit fails, use normal distributions instead
    except RuntimeError:
        popt1, pcov1 = curve_fit(norm.pdf, ofe_arr, odf, p0=(0.1, 0.05))
        fit1 = norm.pdf(ofe_arr, *popt1)
    # reduced chi-square
    chisq1 = reduced_chisqare(odf, fit1, N=popt1.shape[0])
    # KL divergence
    kldiv1 = kl_divergence(odf, fit1, dx)
    
    # two-Gaussian fit
    try:
        popt2, pcov2 = curve_fit(double_skewnorm, ofe_arr, odf, 
                                 p0=(0, 0.1, 0.05, 0.5, 0, 0.3, 0.05),
                                 bounds=([-np.inf, -np.inf, 0., 0., -np.inf, -np.inf, 0.], 
                                         [np.inf, np.inf, np.inf, 1., np.inf, np.inf, np.inf]))
        fit2 = double_skewnorm(ofe_arr, *popt2)
        # normalization coefficient of low-alpha component
        coeff1 = popt2[3]
        coeff1 = max(coeff1, 1 - coeff1)
        # difference in means relative to standard deviation
        delta_mean = np.abs(popt2[5] - popt2[1]) / popt2[2]
    except RuntimeError:
        popt2, pcov2 = curve_fit(double_gaussian, ofe_arr, odf, 
                                 p0=(0.1, 0.05, 0.5, 0.3, 0.05))
        fit2 = double_gaussian(ofe_arr, *popt2)
        # normalization coefficient of low-alpha component
        coeff1 = popt2[2]
        coeff1 = max(coeff1, 1 - coeff1)
        # difference in means relative to standard deviation
        delta_mean = np.abs(popt2[3] - popt2[0]) / popt2[1]
    # reduced chi-square
    chisq2 = reduced_chisqare(odf, fit2, N=popt2.shape[0])
    # KL divergence
    kldiv2 = kl_divergence(odf, fit2, dx)
    # difference in chi-square
    delta_chisq = chisq1 - chisq2
    # difference in KL divergence
    delta_kld = kldiv1 - kldiv2
    
    # diagnostic plot
    if diagnostic:
        fig, ax = plt.subplots()
        ax.plot(ofe_arr, odf, label='VICE')
        ax.plot(ofe_arr, fit1, '--', 
                label='Single peak fit\n(Chi2=%.01e, KLdiv=%.02f)' % (chisq1, kldiv1))
        ax.plot(ofe_arr, fit2, '--', 
                label='Double peak fit\n(Chi2=%.01e, KLdiv=%.02f)' % (chisq2, kldiv2))
        ax.text(0.68, 0.72, 
                r'$\Delta\chi^2=%.02f$' % delta_chisq, 
                transform=ax.transAxes)
        ax.text(0.68, 0.68,
                r'$\Delta D_{\rm{KL}}=%.02f$' % delta_kld,
                transform=ax.transAxes)
        ax.text(0.68, 0.64, 
                r'$\Delta\mu=%.02f\sigma_1$' % delta_mean, 
                transform=ax.transAxes)
        ax.text(0.68, 0.60, 
                r'$C_1=%.02f$' % coeff1, 
                transform=ax.transAxes)
        ax.legend(loc='upper right')
        ax.set_xlabel('[O/Fe]')
        ax.set_ylabel('PDF')
        figname = '_'.join(output_name.split('/')[1:])
        plt.savefig(paths.figures / 'bimodality' / 'skewnorm' / ('%s.png' % figname), dpi=140)
        plt.close()
    
    return delta_chisq, delta_kld, delta_mean, coeff1


def reduced_chisqare(obs, exp, N=1):
    """
    Calculate the reduced chi-square statistic.
    
    Parameters
    ----------
    obs : array-like
        Observed values
    exp : array-like
        Expected values
    N : int, optional
        Number of parameters. The default is 1.
    
    Returns
    -------
    redchisq : float
        Reduced chi-square value.
    """
    # mask zeroes
    obs = obs[exp > 0]
    exp = exp[exp > 0]
    # normalize arrays
    exp *= np.sum(obs) / np.sum(exp)
    chisq = chisquare(obs, exp)[0]
    return chisq / (obs.shape[0] - N)


def double_skewnorm(x, a1, loc1, scale1, coeff1, a2, loc2, scale2):
    """
    A normalized double skew-normal distribution.
    
    Parameters
    ----------
    x : float
        The independent variable.
    a1 : float
        Skewness of the first Gaussian component.
    loc1 : float
        Center of the first Gaussian.
    scale1 : float
        Standard deviation of the first Gaussian.
    coeff1 : float
        Normalization coefficient of the first Gaussian. Also defines the 
        second, which is (1 - coeff1).
    a2 : float
        Skewness of the second Gaussian.
    loc2 : float
        Center of the second Gaussian.
    scale2 : float
        Standard deviation of the second Gaussian.
    """
    func1 = skewnorm.pdf(x, a1, loc1, scale1) 
    func2 = skewnorm.pdf(x, a2, loc2, scale2)
    return coeff1 * func1 + (1 - coeff1) * func2


def double_gaussian(x, loc1, scale1, coeff1, loc2, scale2):
    """
    A normalized double Gaussian distribution.
    
    Parameters
    ----------
    x : float
        The independent variable.
    loc1 : float
        Center of the first Gaussian.
    scale1 : float
        Standard deviation of the first Gaussian.
    coeff1 : float
        Normalization coefficient of the first Gaussian. Also defines the 
        second, which is (1 - coeff1).
    loc2 : float
        Center of the second Gaussian.
    scale2 : float
        Standard deviation of the second Gaussian.
    """
    func1 = norm.pdf(x, loc1, scale1) 
    func2 = norm.pdf(x, loc2, scale2)
    return coeff1 * func1 + (1 - coeff1) * func2


if __name__ == '__main__':
    main()
