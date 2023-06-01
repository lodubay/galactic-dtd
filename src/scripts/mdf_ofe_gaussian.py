"""
This script compares the metallicity distribution functions (MDFs) for 
h277 analog and Gaussian migration schemes.
"""

from ofe_distribution import plot_multiple_comparison

SFH = 'conroy22_JW20yields'
# SFH = 'insideout'
DTD = 'plateau_width300_slope11'
MIGR_LIST = ['diffusion', 'gaussian']
LABELS = ['h277 analog', 'Gaussian migration']

outputs = ['%s/%s/%s' % (mig, SFH, DTD) for mig in MIGR_LIST]
plot_multiple_comparison(outputs, LABELS, fname='mdf_ofe_gaussian.pdf', 
                         verbose=True)
