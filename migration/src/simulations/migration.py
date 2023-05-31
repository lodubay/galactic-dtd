import random
import math as m
from vice.toolkit import hydrodisk
from .._globals import END_TIME, ZONE_WIDTH


class diskmigration(hydrodisk.hydrodiskstars):

    r"""
    A ``hydrodiskstars`` object which writes extra analog star particle data to
    an output file.

    Parameters
    ----------
    radbins : array-like
        The bins in galactocentric radius in kpc corresponding to each annulus.
    mode : ``str`` [default : "linear"]
        A keyword denoting the time-dependence of stellar migration.
        Allowed values:

        - "diffusion"
        - "linear"
        - "sudden"
        - "post-process"

    filename : ``str`` [default : "stars.out"]
        The name of the file to write the extra star particle data to.

    Attributes
    ----------
    write : ``bool`` [default : False]     
        A boolean describing whether or not to write to an output file when
        the object is called. The ``multizone`` object, and by extension the
        ``milkyway`` object, automatically switch this attribute to True at the
        correct time to record extra data.
    """

    def __init__(self, radbins, mode = "diffusion", filename = "stars.out",
        **kwargs):
        super().__init__(radbins, mode = mode, **kwargs)
        if isinstance(filename, str):
            self._file = open(filename, 'w')
            self._file.write("# zone_origin\ttime_origin\tanalog_id\tzfinal\n")
        else:
            raise TypeError("Filename must be a string. Got: %s" % (
                type(filename)))

        # use only disk stars in these simulations
        self.decomp_filter([1, 2])

        # Multizone object automatically swaps this to True in setting up
        # its stellar population zone histories
        self.write = False

    def __call__(self, zone, tform, time):
        if tform == time:
            super().__call__(zone, tform, time) # reset analog star particle
            if self.write:
                if self.analog_index == -1:
                    # finalz = 100
                    finalz = 0
                    analog_id = -1
                else:
                    finalz = self.analog_data["zfinal"][self.analog_index]
                    analog_id = self.analog_data["id"][self.analog_index]
                self._file.write("%d\t%.2f\t%d\t%.2f\n" % (zone, tform,
                    analog_id, finalz))
            else: pass
            return zone
        else:
            return super().__call__(zone, tform, time)

    def close_file(self):
        r"""
        Closes the output file - should be called after the multizone model
        simulation runs.
        """
        self._file.close()

    @property
    def write(self):
        r"""
        Type : bool

        Whether or not to write out to the extra star particle data output
        file. For internal use by the vice.multizone object only.
        """
        return self._write

    @write.setter
    def write(self, value):
        if isinstance(value, bool):
            self._write = value
        else:
            raise TypeError("Must be a boolean. Got: %s" % (type(value)))


class gaussian_migration:
    r"""
    Inherits file read/write functionality from ``diskmigration``.
    """
    def __init__(self, radbins, zone_width=ZONE_WIDTH, end_time=END_TIME,
                 filename="stars.out"):
        self.radial_bins = radbins
        self.zone_width = zone_width
        self.end_time = end_time
        # super().__init__(radbins, mode=None, filename=filename, **kwargs)
        if isinstance(filename, str):
            self._file = open(filename, 'w')
            # Same format as above for compatibility
            self._file.write("# zone_origin\ttime_origin\tanalog_id\tzfinal\n")
        else:
            raise TypeError("Filename must be a string. Got: %s" % (
                type(filename)))

        # Multizone object automatically swaps this to True in setting up
        # its stellar population zone histories
        self.write = False
        
    def __call__(self, zone, tform, time):
        Rform = self.zone_width * (zone + 0.5)
        age = self.end_time - tform
        if tform == time:
            if age > 0:
                # Randomly draw migration distance dR based on age & Rform
                while True: # ensure 0 < Rfinal <= Rmax
                    dR = random.gauss(mu=0., sigma=self.migr_scale(age, Rform))
                    Rfinal = Rform + dR
                    if Rfinal > 0. and Rfinal <= self.radial_bins[-1]:
                        break
            else:
                # Stars born in the last simulation timestep won't migrate
                dR = 0.
                Rfinal = Rform
            self.dR = dR
            # Randomly draw final midplane distance and write to file
            if self.write:
                # list of midplane distances 
                zlist = [i * 0.01 for i in range(-500, 501)]
                # vertical scale height based on age and Rfinal
                hz = self.scale_height(age, Rfinal)
                # draw from sech-squared distribution
                finalz = random.choices(zlist, k=1, 
                    cum_weights=self.sech_squared_cdf(zlist, hz))
                analog_id = -1
                self._file.write("%d\t%.2f\t%d\t%.2f\n" % (zone, tform,
                    analog_id, finalz))
            else:
                pass
            return zone
        else: 
            # Interpolate between Rform and Rfinal at current time
            R = self.interpolator(Rform, Rform + self.dR, tform, time)
            return int((R / self.zone_width) - 0.5)

    def close_file(self):
        r"""
        Closes the output file - should be called after the multizone model
        simulation runs.
        """
        self._file.close()

    @property
    def write(self):
        r"""
        Type : bool

        Whether or not to write out to the extra star particle data output
        file. For internal use by the vice.multizone object only.
        """
        return self._write

    @write.setter
    def write(self, value):
        if isinstance(value, bool):
            self._write = value
        else:
            raise TypeError("Must be a boolean. Got: %s" % (type(value)))
    
    def interpolator(self, Rform, Rfinal, tform, time):
        r"""
        Interpolate between the formation and final radius following the fit
        to the h277 data, $\Delta R \propto t^{0.33}$.
        
        Parameters
        ----------
        Rform : float
            Radius of formation in kpc.
        Rfinal : float
            Final radius in kpc.
        tform : float
            Formation time in Gyr.
        time : float
            Simulation time in Gyr.
        
        Returns
        -------
        float
            Radius at the current simulation time in kpc.
        """
        tfrac = (time - tform) / (self.end_time - tform)
        return Rform + (Rfinal - Rform) * (tfrac ** 0.33)
        
    @staticmethod
    def migr_scale(age, Rform):
        r"""
        A prescription for $\sigma_{\Delta R}$, the scale of the Gaussian 
        distribution of radial migration.
        
        $$ \sigma_{\Delta R} = 1.35 (R_\rm{form}/8\,\rm{kpc})^{0.61} 
        (\tau/1\,\rm{Gyr})^{0.33} $$
        
        Parameters
        ----------
        age : float or array-like
            Age of the stellar population in Gyr.
        Rform : float or array-like
            Formation radius of the stellar population in kpc.
        
        Returns
        -------
        float or array-like
            Scale factor for radial migration $\sigma_{\Delta R}$.
        """
        return 1.35 * (age ** 0.33) * (Rform / 8) ** 0.61
    
    @staticmethod
    def scale_height(age, Rfinal):
        r"""
        The scale height $h_z$ as a function of age and final radius:
        
        $$ h_z = (0.25\,\rm{kpc}) 
        \exp\Big(\frac{(\tau-5\,\rm{Gyr})}{7.0\,\rm{Gyr}}\Big) 
        \exp\Big(\frac{(R_{\rm{final}}-8\,\rm{kpc})}{6.0\,\rm{kpc}}\Big) $$
        
        Parameters
        ----------
        age : float or array-like
            Age in Gyr.
        Rfinal : float or array-like
            Final radius $R_{\rm{final}}$ in kpc.
        
        Returns
        -------
        float
            Scale height $h_z$ in kpc.
        """
        return 0.25 * m.exp((age - 5.) / 7.) * m.exp((Rfinal - 8.) / 6.)
        # return 0.18 * (age ** 0.63) * (Rform / 8) ** 1.15
        
    @staticmethod
    def sech_squared_cdf(x, scale):
        r"""
        The cumulative distribution function (CDF) of the hyperbolic sec-square
        probability distribution function (PDF), which determines the density
        of stars as a function of distance from the midplane $z$. For some
        scaling $h_z$, the PDF is
        
        $$ \rm{PDF}(z) = \frac{1}{4 h_z} \cosh^{-2}\Big(\frac{z}{2 h_z}\Big) $$
        
        and the CDF is
        
        $$ \rm{CDF}(z) = \frac{1}{1 + e^{-z / h_z}} $$
        
        Parameters
        ----------
        x : float
            Independent variable.
        scale : float
            Width of the sech-squared distribution, with the same units as x.
        
        Returns
        -------
        float
            The value of the CDF at the given x.
        """
        return 1 / (1 + m.exp(-x / scale))
        