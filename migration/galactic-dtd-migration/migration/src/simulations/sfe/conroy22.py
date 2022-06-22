r"""
This file declares the time-dependence of the star formation efficiency
timescale from Conroy et al. (2022)
"""

class conroy22:
	r"""
	The piecewise time-dependent star formation efficiency (SFE) timescale model
	from Conroy et al. (2022).

	Parameters
	----------
	t1 : float, optional
		Start time in Gyr of the gradual increase in SFE. The default is 2.5 Gyr.
	t2 : float, optional
		End time in Gyr of the gradual increase in SFE. The default is 3.7 Gyr.
	slope : float, optional
		The steepness of the increase in SFE. The default is 3.
	tau_star_init : float, optional
		The initial SFE timescale in Gyr. The default is 50 Gyr.
	tau_star_final : float, optional
		The final SFE timescale in Gyr. The default is 2.36 Gyr.

	Calling
	-------
	- Parameters

		time : float
			Simulation time in Gyr.
	"""
	def __init__(self, t1 = 2.5, t2 = 3.7, slope = 3, tau_star_init = 50,
			  tau_star_final = 2.36):
		self.t1 = t1
		self.t2 = t2
		self.slope = slope
		self.tau_star_init = tau_star_init
		self.tau_star_final = tau_star_final

	def __call__(self, time, arg2):
		if time < self.t1:
			return self.tau_star_init
		elif time >= self.t1 and time <= self.t2:
			return self.tau_star_init / ((1 + self.slope * (time - self.t1))**2)
		else:
			return self.tau_star_final
