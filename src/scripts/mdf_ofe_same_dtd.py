"""
This script plots multiple [O/Fe] distribution functions (MDFs) for 
multizone runs which assume the same Type Ia delay time distribution (DTD).
"""

from ofe_distribution import plot_multiple_comparison

DTD_NAME = 'powerlaw_slope11'
SFH_LIST = ['insideout',
            'lateburst',
            'conroy22',
            'twoinfall']
LABELS = ['Inside-Out',
          'Late-Burst',
          'Conroy+ 2022',
          'Two-Infall']
MIGRATION = 'diffusion'

def main():
    outputs = ['%s/%s/%s' % (MIGRATION, sfh, DTD_NAME) for sfh in SFH_LIST]
    plot_multiple_comparison(outputs, LABELS, fname='mdf_ofe_same_dtd.pdf', 
                             verbose=True)

if __name__ == '__main__':
    main()
