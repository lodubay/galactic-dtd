import math as m
import numpy as np
import matplotlib.pyplot as plt
import vice
import paths
from multizone.src.simulations import models
from multizone.src.simulations.yields import C22
from track_and_mdf import plot_vice_onezone

area = m.pi * (9**2 - 7**2)
sz = vice.singlezone(name='../data/onezone/conroy22/conroy22',
                     func=models.exponential_ifrmode(8, dr=1), mode='ifr',
                     tau_star=models.conroy22_tau_star(1))
simtime = np.arange(0, 13.21, 0.01)
sz.run(simtime, overwrite=True)

fig, axs = plot_vice_onezone('../data/onezone/conroy22/conroy22')
axs[0].set_xlim((-3, 0.5))
axs[0].set_ylim((-0.1, 0.65))
plt.savefig(paths.figures / 'onezone_conroy22.png', dpi=300)
plt.close()
