"""
This script tests for bimodality in the [O/Fe] distribution of a VICE model
output by fitting a Gaussian to the distribution of an inner and outer region
of the Galaxy and calculating the difference between them.
"""

from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
from utils import multioutput_to_pandas, get_bin_centers
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
        columns=['Significance'],
    )
    
    with tqdm(total=len(SFHs) * len(DTDs)) as t:
        for evolution in SFHs:
            for RIa in DTDs:
                output_name = '/'.join(('diffusion', evolution, RIa))
                summary_table.loc[RIa, evolution] = bimodality(output_name)
                t.update()
    
    print(summary_table)
    summary_table.to_csv(paths.figures / 'bimodality/bimodality.csv')
    

def bimodality(output_name, diagnostic=True):
    """
    Determine if the [O/Fe] distribution of the given VICE output is bimodal.
    
    This function fits a Gaussian to the [O/Fe] distribution in two galactic
    regions: an extreme inner and extreme outer region, both of which should
    be far from the midplane to maximize the fraction of high-alpha stars.
    If the Gaussian fits are substantially different, the distribution is
    bimodal.
    
    Parameters
    ----------
    output_name : str
        Name of VICE multizone output, e.g. 
        'diffusion/insideout/powerlaw_slope11'
    diagnostic : bool
        If True, make a diagnostic plot illustrating the inner and outer ODFs
        and their respective fits.
    
    Returns
    -------
    delta_mean : float
        Difference between the means of the inner and outer region fits
        divided by the RMS standard deviation in data units.
    """
    vice_stars = multioutput_to_pandas(output_name)
    # inner region ODF
    odf_inner, bin_edges = vice_mdf(vice_stars, col='[o/fe]', 
                                    galr_lim=(3, 5), absz_lim=(1, 2), 
                                    xlim=OFE_LIM, bin_width=BIN_WIDTH, 
                                    smoothing=False)
    # outer region ODF
    odf_outer, bin_edges = vice_mdf(vice_stars, col='[o/fe]', 
                                    galr_lim=(13, 15), absz_lim=(0, 0.5), 
                                    xlim=OFE_LIM, bin_width=BIN_WIDTH, 
                                    smoothing=False)
    ofe_arr = get_bin_centers(bin_edges)
    
    # inner region fit
    popt1, pcov1 = curve_fit(norm.pdf, ofe_arr, odf_inner, p0=(0.3, 0.05))
    fit1 = norm.pdf(ofe_arr, *popt1)
    # outer region fit
    popt2, pcov2 = curve_fit(norm.pdf, ofe_arr, odf_outer, p0=(0.1, 0.05))
    fit2 = norm.pdf(ofe_arr, *popt2)
    
    # combine standard deviations in quadrature
    rms_sigma = np.sqrt(0.5 * (popt1[1]**2 + popt2[1]**2))
    center_diff = popt1[0] - popt2[0]
    delta_mean = (center_diff / rms_sigma)
    
    # diagnostic plot
    if diagnostic:
        fig, ax = plt.subplots()
        ax.plot(ofe_arr, odf_inner, label='VICE inner region')
        ax.plot(ofe_arr, odf_outer, label='VICE outer region')
        ax.plot(ofe_arr, fit1, '--', label='Inner region fit')
        ax.plot(ofe_arr, fit2, '--', label='Outer region fit')
        ax.text(0.65, 0.75, 
                r'$\Delta\mu=%.01f\sigma_{\rm{RMS}}$' % delta_mean, 
                transform=ax.transAxes)
        ax.legend(loc='upper right')
        ax.set_xlabel('[O/Fe]')
        ax.set_ylabel('PDF')
        figname = '_'.join(output_name.split('/')[1:])
        plt.savefig(paths.figures / 'bimodality' / ('%s.png' % figname), dpi=140)
        plt.close()
    
    return delta_mean


if __name__ == '__main__':
    main()
