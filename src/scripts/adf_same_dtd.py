"""
This script plots multiple age distribution functions (ADFs) for multizone
runs which assume the same Type Ia delay time distribution (DTD).
"""

from age_distribution import plot_multiple_comparison

DTD_NAME = 'powerlaw_slope11'
SFH_LIST = ['insideout',
            'lateburst',
            'conroy22',
            'twoinfall']
LABELS = ['Inside-Out',
          'Late-Burst',
          'Conroy+ 2022',
          'Two-Infall']

def main(verbose=False):
    outputs = ['diffusion/%s/%s' % (sfh, DTD_NAME) for sfh in SFH_LIST]
    plot_multiple_comparison(outputs, LABELS, fname='adf_same_dtd.png', 
                             verbose=verbose)

if __name__ == '__main__':
    main(verbose=True)
