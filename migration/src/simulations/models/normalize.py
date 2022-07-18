r"""
This file implements the normalization calculation in Appendix A of
Johnson et al. (2021).
"""

from .conroy22_tau_star import conroy22_tau_star
from ..._globals import MAX_SF_RADIUS, END_TIME, M_STAR_MW
import vice
from vice.toolkit import J21_sf_law
import math as m


def normalize(time_dependence, radial_gradient, radius, dt = 0.01, dr = 0.5,
	recycling = 0.4):
	r"""
	Determine the prefactor on the surface density of star formation as a
	function of time as described in Appendix A of Johnson et al. (2021).

	Parameters
	----------
	time_dependence : <function>
		A function accepting time in Gyr and galactocentric radius in kpc, in
		that order, specifying the time-dependence of the star formation
		history at that radius. Return value assumed to be unitless and
		unnormalized.
	radial_gradient : <function>
		A function accepting galactocentric radius in kpc specifying the
		desired stellar radial surface density gradient at the present day.
		Return value assumed to be unitless and unnormalized.
	radius : real number
		The galactocentric radius to evaluate the normalization at.
	dt : real number [default : 0.01]
		The timestep size in Gyr.
	dr : real number [default : 0.5]
		The width of each annulus in kpc.
	recycling : real number [default : 0.4]
		The instantaneous recycling mass fraction for a single stellar
		population. Default is calculated for the Kroupa IMF [1]_.

	Returns
	-------
	A : real number
		The prefactor on the surface density of star formation at that radius
		such that when used in simulation, the correct total stellar mass with
		the specified radial gradient is produced.

	Notes
	-----
	This function automatically adopts the desired maximum radius of star
	formation, end time of the model, and total stellar mass declared in
	``src/_globals.py``.

	.. [1] Kroupa (2001), MNRAS, 322, 231
	"""

	time_integral = 0
	for i in range(int(END_TIME / dt)):
		time_integral += time_dependence(i * dt) * dt * 1.e9 # yr to Gyr

	radial_integral = 0
	for i in range(int(MAX_SF_RADIUS / dr)):
		radial_integral += radial_gradient(dr * (i + 0.5)) * m.pi * (
			(dr * (i + 1))**2 - (dr * i)**2
		)

	return M_STAR_MW / ((1 - recycling) * radial_integral * time_integral)


def normalize_ifrmode(time_dependence, radial_gradient, radius, dt = 0.01,
	dr = 0.1, recycling = 0.4):
	area = m.pi * ((radius + dr)**2 - radius**2)
	tau_star = J21_sf_law(area)
	eta = vice.milkyway.default_mass_loading(radius)
	mgas = 0
	time = 0
	sfh = []
	times = []
	while time < END_TIME:
		sfr = mgas / tau_star(time, mgas) # msun / Gyr
		mgas += time_dependence(time) * dt * 1.e9 # yr-Gyr conversion
		mgas -= sfr * dt * (1 + eta - recycling)
		sfh.append(1.e-9 * sfr)
		times.append(time)
		time += dt
	sfh = vice.toolkit.interpolation.interp_scheme_1d(times, sfh)
	return normalize(sfh, radial_gradient, radius, dt = dt, dr = dr,
		recycling = recycling)


def normalize_conroy22(time_dependence, radial_gradient, radius, dt = 0.01,
	dr = 0.1, recycling = 0.4):
	"""
    A modified version of normalize_ifrmode for use with the Conroy+ 22
    SFE timescale.
    """
	tau_star = conroy22_tau_star()
	eta = vice.milkyway.default_mass_loading(radius)
	mgas = 0
	time = 0
	sfh = []
	times = []
	while time < END_TIME:
		sfr = mgas / tau_star(time, mgas) # msun / Gyr
		mgas += time_dependence(time) * dt * 1.e9 # yr-Gyr conversion
		mgas -= sfr * dt * (1 + eta - recycling)
		sfh.append(1.e-9 * sfr)
		times.append(time)
		time += dt
	sfh = vice.toolkit.interpolation.interp_scheme_1d(times, sfh)
	return normalize(sfh, radial_gradient, radius, dt = dt, dr = dr,
		recycling = recycling)
