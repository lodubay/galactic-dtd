"""
Global variables for the multizone simulations and plotting scripts.
"""

# Simulation globals are located at src/scripts/simulations/_globals.py
from simulations._globals import END_TIME, ZONE_WIDTH, DT

# Plot bins in galactocentric radius and absolute z-height in kpc
GALR_BINS = [3, 5, 7, 9, 11, 13, 15]
ABSZ_BINS = [0, 0.5, 1, 2]

# AASTeX plot widths in inches
ONE_COLUMN_WIDTH = 3.25
TWO_COLUMN_WIDTH = 7.
