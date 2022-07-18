"""
Plot MDFs from all multizone runs
"""

from mdf_9panel import main

for evolution in ['conroy22', 'insideout', 'lateburst']:
    for RIa in ['powerlaw', 'powerlaw_delayed', 'powerlaw_steep',
                'powerlaw_steep_delayed', 'powerlaw_broken', 'bimodal',
                'exponential', 'exponential_delayed', 'exponential_long']:
        print('/'.join((evolution, RIa)))
        main(evolution, RIa)