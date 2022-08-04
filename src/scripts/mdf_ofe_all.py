"""
Plot MDFs from all multizone runs
"""

from mdf_ofe_9panel import main

for evolution in ['insideout_conroy22', 'insideout_johnson21',
                  'lateburst_conroy22', 'lateburst_johnson21',
                  'twoinfall']:
    for RIa in ['powerlaw', 'powerlaw_delayed', 'powerlaw_steep',
                'powerlaw_steep_delayed', 'powerlaw_broken', 'bimodal',
                'exponential', 'exponential_delayed', 'exponential_long']:
        print('/'.join((evolution, RIa)))
        main(evolution, RIa)