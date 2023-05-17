"""
This script plots the change in radius as a function of birth radius for stars
in the h277 hydrodynamical simulation.
"""

import numpy as np
import pandas as pd
from vice.toolkit.hydrodisk import hydrodiskstars
from _globals import ZONE_WIDTH, END_TIME

def main():
    radial_bins = np.arange(0, 20.1, ZONE_WIDTH)
    h277 = hydrodiskstars(radial_bins, N=3102519)
    # Select disk stars only
    h277.decomp_filter([1, 2])
    data = pd.DataFrame(h277.analog_data)
    h277.analog_data["age"] = [END_TIME - _ for _ in h277.analog_data["tform"]]


if __name__ == '__main__':
    main()
