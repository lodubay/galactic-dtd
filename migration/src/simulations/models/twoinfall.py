r"""
This file declares the time-dependence of the star formation history at a
given radius under the two-infall model.
"""

import numpy as np
from numpy.polynomial.polynomial import Polynomial
from ..._globals import END_TIME
from .utils import double_exponential
from .normalize import normalize_ifrmode, twoinfall_ampratio
from .gradient import gradient
import math as m
import os


class twoinfall(double_exponential):

    def __init__(self, radius, dt = 0.01, dr = 0.1, outflows = True):
        spitoni_params = np.genfromtxt('%s/spitoni_twoinfall.dat' % (
            os.path.abspath(os.path.dirname(__file__))))
        super().__init__(onset = 4, ratio = 0.2) # dummy values
        self.first.timescale = self.polyfit(radius, spitoni_params, 3)
        self.second.timescale = self.polyfit(radius, spitoni_params, 5)
        self.onset = self.polyfit(radius, spitoni_params, 9)
        # self.first.timescale = 0.1
        # self.second.timescale = 4
        # self.first.timescale = 0.5 + radius / 15.5
        # self.second.timescale = 2 + 4 * (radius / 15.5)
        self.ratio = self.amp_ratio(radius)
        # self.ratio = twoinfall_ampratio(self, gradient, radius,
        #     onset = self.onset, dt = dt, dr = dr, outflows = outflows)
        prefactor = normalize_ifrmode(self, gradient, radius, dt = dt,
            dr = dr, outflows = outflows)
        self.first.norm *= prefactor
        self.second.norm *= prefactor
        
    def polyfit(self, radius, params, col):
        fit = Polynomial.fit(params[:,0], params[:,col], deg=2, 
                             w=1/params[:,col+1])
        return fit(radius)

    def amp_ratio(self, radius, thin_scale = 2.5, thick_scale = 2.0):
        r"""
        Compute the ratio of amplitudes of the second to the first infall
        episodes.
        Parameters
        ----------
        radius : ``float``
            Galactocentric radius in kpc.
        thin_scale : ``float`` [default : 2.5]
            Scale radius of the thin disk in kpc. Default value chosen from
            Bland-Hawthorn & Gerhard (2016) [1]_.
        thick_scale : ``float`` [default : 2.0]
            Scale radius of the thick disk in kpc. Default value chosen from
            Bland-Hawthorn & Gerhard (2016).
        Returns
        -------
        ratio : ``float``
            The ratio of the amplitudes of the second infall episode to the
            first required to approximately reproduce the balance of thin disk
            and thick disk stars at that radius. See `Notes`_ below.
        Notes
        -----
        The ratio is calculated according to the following:
        .. math:: \left(\frac{\Sigma_t}{\Sigma_T}\right)
            \frac{\tau_1(1 - e^{-T/\tau_1})}{\tau_2(1 - e^{-(T - t_2)/\tau_2})}
        where :math:`\Sigma_t/Sigma_T` is the ratio of thin-disk to thick-disk
        stars at that radius, :math:`\tau_1` is the e-folding timescale of the
        first accretion event, :math:`\tau_2` is the e-folding timescale of the
        second, :math:`t_2` is the time of onset of the second episode, and
        :math:`T` is the total time the simulation runs for. The ratio of
        thin-disk to thick-disk stars is compute according to the following:
        .. math:: \frac{1}{0.27} exp\left(R\frac{R_t - R_T}{R_t R_T}\right)
        where :math:`R` is the Galactocentric radius, :math:`R_t` is the
        scale radius of the thin disk, :math:`R_T` is the scale radius of the
        thick disc, and :math:`1 / 0.27` is their ratio at :math:`R = 0`. This
        can be derived algebraically from relating the integral over the two
        exponential components of the infall history to the ratio of thin-disk
        to thick-disk stars.
        .. [1] Bland-Hawthorn & Gerhard (2016), ARA&A, 54, 529
        """
        thin_to_thick = m.exp(radius * (1 / thick_scale - 1 / thin_scale))
        # thin_to_thick = m.exp(radius * (thin_scale - thick_scale) / (
        #     thin_scale * thick_scale))
        thin_to_thick /= 0.27
        timescale_factor = self.first.timescale / self.second.timescale
        timescale_factor *= (1 - m.exp(-END_TIME / self.first.timescale))
        timescale_factor /= (1 - m.exp(-(END_TIME - self.onset) / 
                                       self.second.timescale))
        return thin_to_thick * timescale_factor
