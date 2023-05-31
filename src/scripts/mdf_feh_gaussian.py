"""
This script compares the metallicity distribution functions (MDFs) for 
h277 analog and Gaussian migration schemes.
"""

from feh_distribution import plot_multiple_comparison

SFH = 'insideout'
DTD = 'exponential_timescale15'
MIGR_LIST = ['diffusion', 'gaussian']
LABELS = ['h277 analog', 'Gaussian migration']

outputs = ['%s/%s/%s' % (mig, SFH, DTD) for mig in MIGR_LIST]
plot_multiple_comparison(outputs, LABELS, fname='mdf_feh_gaussian.pdf', 
                         verbose=True)
