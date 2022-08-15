"""
This script plots multiple age distribution functions (ADFs) for multizone
runs which assume the same star formathion history (SFH).
"""

from age_distribution import plot_multiple_comparison

SFH_NAME = 'insideout_johnson21'
DTD_LIST = ['powerlaw_slope11_delay040',
            'exponential_timescale30',
            'plateau_width200_slope11',
            'prompt_peak050_timescale30']
LABELS = ['Power-Law\n($\\alpha=-1.1$)',
          'Exponential\n($\\tau=3$ Gyr)',
          'Broken Power-Law\n($W=200$ Myr)',
          'Prompt + Exponential\n($\\mu=50$ Myr)']

def main(verbose=False):
    outputs = ['diffusion/%s/%s' % (SFH_NAME, dtd) for dtd in DTD_LIST]
    plot_multiple_comparison(outputs, LABELS, fname='adf_same_sfh.pdf', 
                             verbose=verbose, double_line_titles=True)

if __name__ == '__main__':
    main(verbose=True)
