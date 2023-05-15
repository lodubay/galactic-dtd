"""
This script calculates a numerical score for how well a given VICE output
matches the APOGEE data in the Age-[O/Fe] plane. For a single [O/Fe] bin,
the difference in the medians is divided by the RMS 1-sigma bounds of the
distribution. A given Galactic region is scored by the RMS of the difference
between medians, weighted by the number of APOGEE stars per bin. The overall
score for the model is the average of the scores for each region weighted
by the number of APOGEE stars in each region.
"""

from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scatter_plot_grid import setup_axes
from utils import import_apogee, apogee_region, multioutput_to_pandas, \
    filter_multioutput_stars, weighted_quantile, group_by_bins, \
    get_bin_centers, model_uncertainties
from _globals import GALR_BINS, ABSZ_BINS, ZONE_WIDTH
import paths

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

GALR_BINS = GALR_BINS[:-1]
# AGE_LIM_LINEAR = (-1, 14)
# AGE_LIM_LOG = (0.2, 20)
OFE_LIM = (-0.15, 0.55)
OFE_BIN_WIDTH = 0.05
AGE_SOURCES = ['F18', # Feuillet et al. 2018, solar neighborhood only
               'M19', # Mackereth et al. 2019, astroNN
               'L23'] # Leung et al. 2023, variational encoder-decoder
AGE_LABELS = {'F18': 'Feuillet et al. 2018',
              'M19': 'Mackereth et al. 2019',
              'L23': 'Leung et al. 2023'}

def main():
    """
    Score all models provided in the lists above and export a summary table.
    """
    apogee_data = import_apogee(verbose=True)
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTDs, SFHs], names=['DTD', 'SFH']),
        columns=['Score'],
    )
    
    with tqdm(total=len(SFHs) * len(DTDs)) as t:
        for evolution in SFHs:
            for RIa in DTDs:
                output_name = '/'.join(('diffusion', evolution, RIa))
                scores = score_model(output_name, apogee_data)
                summary_table.loc[RIa, evolution] = scores
                t.update()
    
    print(summary_table)
    summary_table.to_csv(paths.figures / 'age_ofe/scores/scores.csv')
    

def score_model(output_name, apogee_data):
    """
    Calculate the score for a single VICE model representing how well it
    matches APOGEE data in the Age-[O/Fe] plane.
    """
    vice_stars = multioutput_to_pandas(output_name)
    vice_stars = model_uncertainties(vice_stars.copy())
    
    fig, axs = setup_axes(xlim=(0, 10), ylim=OFE_LIM, 
                          xlabel='Age difference [Gyr]', ylabel='[O/Fe]')
    
    scores = []
    weights = []
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            vice_subset = filter_multioutput_stars(vice_stars, galr_lim, absz_lim)
            apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
            apogee_subset.dropna(axis=0, subset=['O_FE', 'LATENT_AGE'], inplace=True)
            d_rms = rms_median_diff(vice_subset, apogee_subset, ax,
                                    ofe_lim=OFE_LIM, ofe_bin_width=OFE_BIN_WIDTH)
            scores.append(d_rms)
            weights.append(apogee_subset.shape[0])
            del apogee_subset
            del vice_subset
    
    axs[0,-1].legend()
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    figname = '_'.join(output_name.split('/')[1:])
    plt.savefig(paths.figures / ('age_ofe/scores/%s.png' % figname), dpi=300)
    plt.close()
    
    del vice_stars    
    return np.average(scores, weights=weights)
    
    
def rms_median_diff(vice_stars, apogee_data, ax, ofe_lim=(-0.15, 0.55), 
                    ofe_bin_width=0.05):
    """
    Calculate the RMS median difference between VICE and APOGEE ages.
    """
    ofe_bins = np.arange(ofe_lim[0], ofe_lim[1]+ofe_bin_width, ofe_bin_width)
    # bin APOGEE ages by [O/Fe]
    apogee_grouped = group_by_bins(apogee_data, 'O_FE', ofe_bins)['LATENT_AGE']
    apogee_medians = apogee_grouped.median()
    # count all APOGEE stars in each bin
    apogee_counts = apogee_grouped.count()
    # bin mass-weighted VICE ages by [O/Fe]
    vice_grouped = group_by_bins(vice_stars, '[o/fe]', bins=ofe_bins)
    # weighted medians of VICE
    wm = lambda x: weighted_quantile(x, 'age', 'mass', quantile=0.5)
    vice_medians = vice_grouped.apply(wm)
    # RMS of median difference
    notna = (pd.notna(apogee_medians) & pd.notna(vice_medians))
    median_diffs = vice_medians[notna] - apogee_medians[notna]
    d_rms = np.sqrt(np.average(median_diffs**2, weights=apogee_counts[notna]))
    # plot
    ax.scatter(np.abs(median_diffs), get_bin_centers(ofe_bins)[notna],
               s=np.sqrt(apogee_counts[notna] / 10), c='k')
    ax.axvline(d_rms, color='gray', label='Weighted RMS diff')
    
    return d_rms


if __name__ == '__main__':
    main()
