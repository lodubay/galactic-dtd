"""
This script plots the [O/Fe] bimodality for all VICE multizone runs.
"""

from tqdm import tqdm
from ofe_bimodality import main

SFH_LIST = [
    'insideout', 
    'lateburst', 
    'earlyburst', 
    'twoinfall'
]
DTD_LIST = [
    'powerlaw_slope11', 
    'powerlaw_slope14', 
    'exponential_timescale15',
    'exponential_timescale30', 
    'plateau_width03', 
    'plateau_width10', 
    'prompt',
    'triple'
]

with tqdm(total=len(SFH_LIST) * len(DTD_LIST) + 1) as t:
    for evolution in SFH_LIST:
        for RIa in DTD_LIST:
            # Import VICE multi-zone output data
            output_name = '/'.join(['gaussian', evolution, RIa, 'diskmodel'])
            main(output_name, uncertainties=True)
            t.update()
    output_name = 'diffusion/insideout/powerlaw_slope11/diskmodel'
    main(output_name, uncertainties=True)
    t.update()
