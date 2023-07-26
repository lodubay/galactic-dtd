"""
This script generates three plots which demonstrate the effect of changing
DTD parameters on one-zone model outputs.
"""

from onezone_powerlaw_slope import main as plot_powerlaw_slope
from onezone_exponential_timescale import main as plot_exponential_timescale
from onezone_plateau_width import main as plot_plateau_width

plot_powerlaw_slope()
plot_exponential_timescale()
plot_plateau_width()
