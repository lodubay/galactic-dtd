"""
This script generates [O/Fe]-[Fe/H] grid plots for all VICE multi-zone runs.
"""

from ofe_feh_grid import main
from tqdm import tqdm

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
            output_name = '/'.join(('gaussian', evolution, RIa, 'diskmodel'))
            main(output_name, uncertainties=True, tracks=True, contours=True)
            t.update()
    output_name = 'diffusion/insideout/powerlaw_slope11/diskmodel'
    main(output_name, uncertainties=True, tracks=True, contours=True)
    t.update()
