# -*- coding: utf-8 -*-

"""
Setup classes for defining default parameters or loading in parameters from files:
> Spectrometer attributes
"""
import numpy as np


class SpecSpecs:
    """Object containing information on spectrometer setup and acquisition settings

    Parameters
    ----------
    filename : str
        path to configuration file (*.txt) to read in camera parameters. If no file is provided the internal defaults
        are use
    """
    def __init__(self, filename=None):
        self.filename = filename    # Filename for loading specifications

        # Hidden variables
        self._bit_depth = None  # Hidden bit depth holder
        self._max_DN = None     # Maximum digital number of images (associated with bit depth)
        self._int_time_idx = None

        self._default_specs()       # Load default specs to start

        # Load configuration from text file if it is provided
        if isinstance(self.filename, str):
            self.load_specs(self.filename)

    def _default_specs(self):
        """Define spectrometer default specs > Flame-S"""
        # Spectrometer specs
        self.model = "Flame-S"      # Spectrometer model
        self.fov = None             # Field of view fo spectrometer
        self.ILS = None             # Number array holding instrument line shape (possibly don't hold this here?)
        self.pix_num = 2048         # Number of pixels
        self.bit_depth = 16         # Bit depth of spectrometer detector

        # File information
        self.file_ext = '.npy'  # Spectra saved as numpy array
        self.file_ss = '{}ss'  # Shutter speed format spec
        self.file_spec_type = {'meas': 'Plume', 'dark': 'Dark', 'cal': 'ppmm', 'clear': 'Clear'}
        self.file_datestr = "%Y-%m-%dT%H%M%S"                   # Date/time format spec in filename


        # Acquisition settings
        self.start_int_time = 100       # Starting integration time
        self.start_coadd = 5            # Number of spectra to coadd
        self.framerate = 1              # Framerate of acquisitions (Hz)
        self.wavelengths = None         # Wavelengths (nm)
        self.spectrum = None            # Spectrum
        self.spectrum_filename = None   # Filename for spectrum

        self.auto_int = True        # Bool for requesting automated integration time adjustment
        self.min_saturation = 0.5   # Minimum saturation accepted before adjusting shutter speed (if auto_ss is True)
        self.max_saturation = 0.9   # Maximum saturation accepted before adjusting shutter speed (if auto_ss is True)
        self.saturation_range = [320, 330]  # Range of wavelengths used in checking integration time

        # Predefined list of integration times for automatic exposure adjustment
        self.int_list = np.concatenate((np.arange(1, 10, 1),
                                        np.arange(10, 50, 5),
                                        np.arange(50, 100, 10),
                                        np.arange(100, 500, 50),
                                        np.arange(500, 1000, 100),
                                        np.arange(10 ** 3, 10 ** 4, 500),
                                        np.array([10 ** 4])))

    def load_specs(self, filename):
        """Load spectrometer specifications from file

        Parameters
        ----------
        filename : str
            path to configuration (*.txt) file
        """
        self.filename = filename
        # Add loading functionality here

    def save_specs(self, filename):
        """Save spectrometer specifications to file

        Parameters
        ----------
        filename : str
            path to configuration (*.txt) file
        """
        pass

    @property
    def bit_depth(self):
        return self._bit_depth

    @bit_depth.setter
    def bit_depth(self, value):
        """Update _max_DN when bit_depth is defined (two are intrinsically linked)"""
        self._bit_depth = value
        self._max_DN = (2 ** self.bit_depth) - 1
