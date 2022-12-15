"""
Plot #3 for my AAS 241. This plot compares the [O/Fe] distributions from
APOGEE, the inside-out + power-law model, the Conroy22 + power-law model, 
and the Conroy22 + short plateau model.
"""

import matplotlib.pyplot as plt
from ofe_distribution import plot_multiple_comparison

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

outputs = ['diffusion/conroy22/powerlaw_slope11',
           'diffusion/conroy22/plateau_width300_slope11']

labels = ['Power Law',
          'Power Law with Plateau']

plot_multiple_comparison(outputs, labels, verbose=True, fname='presentation_plot03.png')