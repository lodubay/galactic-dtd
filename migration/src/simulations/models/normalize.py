r"""
This file implements the normalization calculation in Appendix A of
Johnson et al. (2021).
"""

from .conroy22_tau_star import conroy22_tau_star
from .twoinfall_tau_star import twoinfall_tau_star
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
                      dr = 0.1, recycling = 0.4, which_tau_star='johnson21',
                      outflows = True):
    r"""
    Performs essentially the same thing as ``normalize`` but for models ran in
    infall mode.
    """
    area = m.pi * ((radius + dr)**2 - radius**2)
    if which_tau_star.lower() == 'spitoni21':
        # tau_star = twoinfall_tau_star(area, radius)
        tau_star = lambda t, mgas: 2.
    elif which_tau_star.lower() == 'conroy22':
        tau_star = conroy22_tau_star(area)
    elif which_tau_star.lower() == 'johnson21':
        tau_star = J21_sf_law(area)
    else:
        raise TypeError('Unrecognized prescription for tau_star!')
    # tau_star = {
    #     'johnson21': J21_sf_law,
    #     'conroy22': conroy22_tau_star,
    #     'spitoni21': spitoni21_tau_star,
    # }[which_tau_star](area)
    if outflows:
        eta = vice.milkyway.default_mass_loading(radius)
    else:
        eta = 0
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


def twoinfall_ampratio(time_dependence, radial_gradient, radius, onset = 4,
                       dt = 0.01, dr = 0.1, recycling = 0.4, thin_scale = 2.0, 
                       thick_scale = 2.5):
    area = m.pi * ((radius + dr)**2 - radius**2)
    tau_star = J21_sf_law(area)
    eta = vice.milkyway.default_mass_loading(radius)
    mgas = 0
    time = 0
    mstar = 0
    mstar_at_onset = None
    while time < END_TIME:
        sfr = mgas / tau_star(time, mgas) # msun / Gyr
        mgas += time_dependence(time) * dt * 1.e9 # yr-Gyr conversion
        mgas -= sfr * dt * (1 + eta - recycling)
        mstar += sfr * dt * (1 - recycling)
        if mstar_at_onset is None and time >= onset: mstar_at_onset = mstar
        time += dt
    thick_to_thin = 0.27 * m.exp(radius * (1 / thin_scale - 1 / thick_scale))
    return mstar / (mstar - mstar_at_onset) * (1 + thick_to_thin)**-1
