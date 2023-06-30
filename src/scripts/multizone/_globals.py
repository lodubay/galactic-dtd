"""
Global variables for the multizone simulations.
"""

# Total simulation time in Gyr
END_TIME = 13.2

# Width of each annulus in kpc
ZONE_WIDTH = 0.1

# Radius in kpc beyond which the SFR = 0
MAX_SF_RADIUS = 15.5

# Stellar mass of Milky Way (Licquia & Newman 2015, ApJ, 806, 96)
M_STAR_MW = 5.17e10

# Simulation timestep size in Gyr
DT = 0.01

# List of elements to simulate
ELEMENTS = ['o', 'fe']

# Number of stellar populations per zone per timestep
NSTARS = 8

# Whether to force overwrite existing VICE outputs of the same name
OVERWRITE = True

# The minimum SN Ia delay time in Gyr
MIN_RIA_DELAY = 0.04
