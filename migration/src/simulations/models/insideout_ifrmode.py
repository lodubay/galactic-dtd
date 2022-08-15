r"""
This file declares the time-dependence of the star formation history at a
given radius in a modified version of the inside-out model from Johnson et al. 
(2021).
"""

from .insideout import insideout
from .normalize import normalize_ifrmode
from .gradient import gradient

class insideout_ifrmode(insideout):

    r"""
    The inside-out SFH model from Johnson et al. (2021) modified to run in
    infall mode.

    Parameters
    ----------
    radius : float
        The galactocentric radius in kpc of a given annulus in the model.
    dt : float [default : 0.01]
        The timestep size of the model in Gyr.
    dr : float [default : 0.1]
        The width of the annulus in kpc.
    which_tau_star : str [default: "johnson21"]
        Which SFE timescale prescription to assume.
        Allowed values:
            
        - "johnson21"
        - "conroy22"

    Functions
    ---------
    - timescale [staticmethod]

    Other atributes and functionality are inherited from
    ``insideout`` declared in ``src/simulations/models/insideout.py``.
    """

    def __init__(self, radius, dt = 0.01, dr = 0.1, which_tau_star = 'johnson21'):
        super().__init__(radius, dt = dt, dr = dr)
        self.norm = normalize_ifrmode(self, gradient, radius, dt = dt, dr = dr,
                                      which_tau_star = which_tau_star)
