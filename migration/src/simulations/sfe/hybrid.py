from vice.toolkit import J21_sf_law
from .conroy22 import conroy22

class hybrid(J21_sf_law, conroy22):

	def __init__(self, area, mode = "sfr", t1 = 2.5, t2 = 3.7, slope = 3,
		tau_star_init = 50, tau_star_final = 2.36, **kwargs):

		J21_sf_law.__init__(self, area, mode = mode, **kwargs)
		conroy22.__init__(self, t1 = t1, t2 = t2, slope = slope,
			tau_star_init = tau_star_init, tau_star_final = tau_star_final)

	def __call__(self, time, arg2):
		return J21_sf_law.__call__(self, time, arg2) * conroy22.__call__(self, time)