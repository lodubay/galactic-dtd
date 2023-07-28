"""
This script compares the results of multi-zone simulations to APOGEE data.
Various statistical tests are used to assign a score for each diagnostic plot:
[Fe/H] distribution functions, [O/Fe] distribution functions, the [Fe/H] -
[O/Fe] plane, the age - [O/Fe] plane, and bimodality of the [O/Fe] distribution.
"""

from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from apogee_tools import import_apogee, apogee_region, apogee_mdf
from multizone_stars import MultizoneStars
from utils import kl_divergence, kl_div_2D, group_by_bins, weighted_quantile
import paths
from _globals import ABSZ_BINS, GALR_BINS

MIGRATION = 'gaussian'
SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
DTD_LIST = ['powerlaw_slope11', 
            'powerlaw_slope14', 
            'exponential_timescale15',
            'exponential_timescale30',
            'plateau_width03',
            'plateau_width10',
            'prompt',
            'triple']
AGE_SOURCE = 'L23' # Use Leung et al (2023) "latent-space" ages

def main():    
    apogee_data = import_apogee()
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTD_LIST, SFH_LIST], 
                                         names=['DTD', 'SFH']),
        # temporary column names
        columns=['feh_df', 'ofe_df', 'ofe_feh', 'age_ofe', 'bimodality'],
    )
    
    age_col = {'L23': 'LATENT_AGE', 'M19': 'ASTRONN_AGE'}[AGE_SOURCE]
    
    with tqdm(total=len(SFH_LIST) * len(DTD_LIST)) as t:
        for dtd in DTD_LIST:
            for sfh in SFH_LIST:
                output_name = '/'.join([MIGRATION, sfh, dtd, 'diskmodel'])
                mzs = MultizoneStars.from_output(output_name)
                mzs.model_uncertainty(apogee_data, inplace=True)
                # Set up lists to track scores by region
                scores = {col: [] for col in summary_table.columns[:-1]}
                weights = []
                for i in range(len(ABSZ_BINS) - 1):
                    absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
                    for j in range(len(GALR_BINS) - 1):
                        galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
                        # Select stars by galactic region
                        vice_subset = mzs.region(galr_lim, absz_lim)
                        apogee_subset = apogee_region(apogee_data, 
                                                      galr_lim, absz_lim)
                        # Calculate scores for different parameter spaces
                        scores['feh_df'].append(
                            score_feh_df(vice_subset, apogee_subset))
                        scores['ofe_df'].append(
                            score_ofe_df(vice_subset, apogee_subset))
                        scores['ofe_feh'].append(
                            score_ofe_feh(vice_subset, apogee_subset))
                        scores['age_ofe'].append(
                            score_age_ofe(vice_subset, apogee_subset, 
                                          age_col=age_col))
                        weights.append(apogee_subset.shape[0])
                # Append weighted mean scores
                weighted_sums = {col: np.average(scores[col], weights=weights) 
                                 for col in summary_table.columns[:-1]}
                # Separate bimodality test (binary output)
                weighted_sums['bimodality'] = test_bimodality(mzs)
                summary_table.loc[dtd, sfh] = weighted_sums
                t.update()
    
    summary_table.to_csv(paths.data / 'multizone' / 'scores.csv')
    

def score_feh_df(mzs, apogee_data, data_range=(-3., 1.), bin_width=0.01):
    """
    Calculate the KL divergence between VICE and APOGEE MDFs in [Fe/H].
    
    Parameters
    ----------
    mzs : MultizoneStars object
        Star particle output from a VICE multizone run.
    apogee_data : pandas.DataFrame
        APOGEE data
    data_range : tuple, optional
        Lower and upper limits on [Fe/H]. The default is (-3., 1.).
    bin_width : float, optional
        Width of histogram bins in dex. The default is 0.01.
    
    Returns
    -------
    float
        The KL divergence, which represents the difference between the MDFs.
        A lower value represents a closer match.
    """
    bin_edges = np.arange(data_range[0], data_range[1] + bin_width, bin_width)
    vice_dist, _ = mzs.mdf('[fe/h]', bins=bin_edges)
    apogee_dist, _ = apogee_mdf(apogee_data, col='FE_H', bins=bin_edges)
    return kl_divergence(apogee_dist, vice_dist, bin_width)


def score_ofe_df(mzs, apogee_data, data_range=(-0.2, 0.6), bin_width=0.005):
    """
    Calculate the KL divergence between VICE and APOGEE MDFs in [O/Fe].
    
    Parameters
    ----------
    mzs : MultizoneStars object
        Star particle output from a VICE multizone run.
    apogee_data : pandas.DataFrame
        APOGEE data
    data_range : tuple, optional
        Lower and upper limits on [O/Fe]. The default is (-0.2, 0.6).
    bin_width : float, optional
        Width of histogram bins in dex. The default is 0.005.
    
    Returns
    -------
    float
        The KL divergence, which represents the difference between the MDFs.
        A lower value represents a closer match.
    """
    bin_edges = np.arange(data_range[0], data_range[1] + bin_width, bin_width)
    vice_dist, _ = mzs.mdf('[o/fe]', bins=bin_edges)
    apogee_dist, _ = apogee_mdf(apogee_data, col='O_FE', bins=bin_edges)
    return kl_divergence(apogee_dist, vice_dist, bin_width)


def score_ofe_feh(mzs, apogee_data):
    """
    Calculate the 2D KL divergence between VICE and APOGEE in the [O/Fe]-[Fe/H]
    parameter space.
    
    Parameters
    ----------
    mzs : MultizoneStars object
        Star particle output from a VICE multizone run.
    apogee_data : pandas.DataFrame
        APOGEE data
        
    Returns
    -------
    float
        The 2D KL divergence, with a lower value representing a closer match
        between the distributions.
    """
    # 2D KL divergence between APOGEE (true) and VICE (approximate)
    return kl_div_2D(apogee_data[['FE_H', 'O_FE']], mzs(['[fe/h]', '[o/fe]']))


def score_age_ofe(mzs, apogee_data, age_col='LATENT_AGE', 
                  ofe_range=(-0.15, 0.55), bin_width=0.05):
    """
    Calculate the RMS of the difference in medians between VICE and data ages.
    
    Parameters
    ----------
    mzs : MultizoneStars object
        Star particle output from a VICE multizone run.
    apogee_data : pandas.DataFrame
        APOGEE data with ages, typically a subset of a galactic region.
    age_col : str, optional
        Name of column with age data in apogee_data. The default is 'LATENT_AGE'
        which is the Leung et al. (2023) ages.
    ofe_range : tuple, optional
        Outermost bounds on [O/Fe]. The default is (-0.15, 0.55).
    bin_width : float, optional
        The [O/Fe] bin width in dex. The default is 0.05.
        
    Returns
    -------
    float
        Root-mean-square of the difference in median ages in each [O/Fe] bin.
    """
    # Error handling
    age_col_options = ['LATENT_AGE', 'ASTRONN_AGE']
    if age_col not in age_col_options:
        raise ValueError('Parameter "age_col" must be one of', age_col_options)
    ofe_bins = np.arange(ofe_range[0], ofe_range[1] + bin_width, bin_width)
    # bin APOGEE ages by [O/Fe]
    apogee_grouped = group_by_bins(apogee_data, 'O_FE', ofe_bins)[age_col]
    apogee_medians = apogee_grouped.median()
    # count all APOGEE stars in each bin
    apogee_counts = apogee_grouped.count()
    # bin mass-weighted VICE ages by [O/Fe]
    vice_grouped = group_by_bins(mzs.stars, '[o/fe]', bins=ofe_bins)
    # weighted medians of VICE
    wm = lambda x: weighted_quantile(x, 'age', 'mass', quantile=0.5)
    vice_medians = vice_grouped.apply(wm)
    # RMS of median difference
    notna = (pd.notna(apogee_medians) & pd.notna(vice_medians))
    median_diffs = vice_medians[notna] - apogee_medians[notna]
    return np.sqrt(np.average(median_diffs**2, weights=apogee_counts[notna]))


def test_bimodality(mzs, prominence=0.1, feh_bins=[(-0.6, -0.4), (-0.4, -0.2)],
                    galr_lim=(7, 9), absz_lim=(0, 2), smoothing=0.05):
    """
    Determine whether the distribution of stars in [O/Fe] is bimodal.
    
    Parameters
    ----------
    mzs : MultizoneStars object
        Star particle output from a VICE multizone run.
    prominence : float, optional
        Prominence threshold for peak-finding algorithm. The default is 0.1.
    feh_bins : list of tuples, optional
        Limits on [Fe/H] to slice the output. Bimodality will be checked for
        stars within each slice, and if either slice is bimodal, the function
        will return True. The default is [(-0.6, -0.4), (-0.4, -0.2)].
    galr_lim : tuple, optional
        Limits on galactic radius in kpc. The default is (7, 9).
    absz_lim : tuple, optional
        Limits on absolute galactic z-height in kpc. The default is (0, 2).
    smoothing : float, optional
        Boxcar smoothing width for the [O/Fe] distribution. The default is 0.05.
    
    Returns
    -------
    bool
        Whether or not the distribution of stars in [O/Fe] is bimodal.
    """
    subset = mzs.region(galr_lim=galr_lim, absz_lim=absz_lim)
    for feh_bin in feh_bins:
        subset_slice = subset.filter({'[fe/h]': feh_bin})
        mdf, bin_edges = subset_slice.mdf('[o/fe]', smoothing=smoothing,
                                          bins=np.arange(-0.15, 0.56, 0.01))
        peaks, _ = find_peaks(mdf/mdf.max(), prominence=prominence)
        is_bimodal = (len(peaks) > 1)
        if is_bimodal:
            break
    return is_bimodal


if __name__ == '__main__':
    main()
