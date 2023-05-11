"""
Global variables for use by multiple plots.
"""

GALR_BINS = [3, 5, 7, 9, 11, 13, 15] # kpc
ABSZ_BINS = [0, 0.5, 1, 2] # kpc
ZONE_WIDTH = 0.1 # kpc
DT = 0.01 # Gyr
END_TIME = 13.2 # Gyr

# AASTeX plot widths
ONE_COLUMN_WIDTH = 3.25
TWO_COLUMN_WIDTH = 7.

# Error fit polynomial degree for APOGEE parameters
ERROR_FIT_DEG = {
    'FE_H': 1,
    'O_FE': 2,
    'LOG_LATENT_AGE': 3
}

# Error fit range for APOGEE parameters
ERROR_FIT_RANGE = {
    'FE_H': None,
    'O_FE': (-0.2, 0.7),
    'LOG_LATENT_AGE': None
}
