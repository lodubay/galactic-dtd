import random
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


class gaussian_migration(diskmigration):
    r"""
    Inherits file read/write functionality from ``diskmigration``.
    """
    def __init__(self, radbins, zone_width=ZONE_WIDTH, filename="stars.out",
                 **kwargs):
        self.zone_width = zone_width
        super().__init__(radbins, mode=None, filename=filename, **kwargs)
        
    def __call__(self, zone, tform, time):
        Rform = self.zone_width * (zone + 0.5)
        age = END_TIME - tform
        if tform == time:
            # Randomly draw migration distance dR based on age & Rform
            while True: # ensure Rfinal > 0
                dR = random.gauss(loc=0, scale=self.migr_scale(age, Rform))
                Rfinal = Rform + dR
                if Rfinal > 0:
                    break
            self.dR = dR
            # Randomly draw final midplane distance and write to file
            if self.write:
                finalz = random.expovariate(1 / self.scale_height(age, Rform))
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
    
    def interpolator(Rform, Rfinal, tform, time):
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
        tfrac = (time - tform) / (END_TIME - tform)
        return Rform + (Rfinal - Rform) * (tfrac ** 0.33)
        
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
    
    def scale_height(age, Rform):
        r"""
        The scale height $h_z$ as a function of age and formation radius:
        
        $$ h_z = 0.18 (R_\rm{form}/8\,\rm{kpc})^{1.15} 
        (\tau/1\,\rm{Gyr})^{0.63} $$
        
        Parameters
        ----------
        age : float or array-like
            Age in Gyr.
        rform : float or array-like
            Formation radius $R_{\rm{form}}$ in kpc.
        
        Returns
        -------
        float
            Scale height $h_z$ in kpc.
        """
        return 0.18 * (age ** 0.63) * (Rform / 8) ** 1.15
        