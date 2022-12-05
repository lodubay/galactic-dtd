"""
Generate a grid of age vs [O/Fe] plots for multiple Galactic regions.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
import vice
from utils import multioutput_to_pandas, filter_multioutput_stars, \
    sample_dataframe
from _globals import GALR_BINS, ABSZ_BINS
from ofe_feh_vice import setup_colorbar
import paths


def main(output_name, data_dir='../data/migration', cmap='winter'):
    pass


if __name__ == "__main__":
    evolution = sys.argv[1]
    RIa = sys.argv[2]
    main('/'.join(['diffusion', evolution, RIa]))
