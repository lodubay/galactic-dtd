"""
This script generates two plots to be included in a single figure: the effect
of the minimum SN Ia delay time and the different DTD models on the one-zone
model outputs.
"""

from onezone_minimum_delay import main as plot_minimum_delay
from onezone_dtd import main as plot_all

plot_minimum_delay()
plot_all()
