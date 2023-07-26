"""
This script generates plots for a multi-panel figure comparing the results
of several one-zone chemical evolution models with various DTD parameters.
"""

from onezone_powerlaw_slope import main as plot_powerlaw_slope
from onezone_exponential_timescale import main as plot_exponential_timescale
from onezone_plateau_width import main as plot_plateau_width
from onezone_minimum_delay import main as plot_minimum_delay
from onezone_dtd import main as plot_all

plot_powerlaw_slope()
plot_exponential_timescale()
plot_plateau_width()
plot_minimum_delay()
plot_all()
