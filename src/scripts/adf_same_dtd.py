"""
This script plots multiple age distribution functions (ADFs) for multizone
runs which assume the same star formathion history (SFH).
"""

from age_distribution import plot_multiple_comparison

DTD_NAME = 'powerlaw_slope11_delay040'
SFH_LIST = ['insideout_johnson21',
            'lateburst_johnson21',
            'insideout_conroy22',
            'twoinfall']
LABELS = ['Inside-Out',
          'Late-Burst',
          'Conroy+ 2022',
          'Two-Infall']

def main(verbose=False):
    outputs = ['diffusion/%s/%s' % (sfh, DTD_NAME) for sfh in SFH_LIST]
    plot_multiple_comparison(outputs, LABELS, fname='adf_same_dtd.pdf', 
                             verbose=verbose)

if __name__ == '__main__':
    main(verbose=True)
