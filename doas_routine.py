# Main Subroutine which processes images according to the DOAS retrieval method.

from automation.setupclasses import SpecSpecs
from automation.utils import load_spectrum
from automation.controllers import ScanProperties
from automation.directory_watcher import create_dir_watcher
from save_spec import SaveSpectra

import numpy as np
from scipy import signal
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from astropy.convolution import convolve
import scipy.integrate as integrate
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
import threading
import queue
from tkinter import filedialog
from tkinter import messagebox
import datetime
import copy

warnings.filterwarnings("ignore", category=OptimizeWarning)


class DOASWorker:
    """
    Class to control DOAS processing
    General order of play for processing:
    Initiate class,
    get_ref_spectrum()
    set_fit_window()
    shift_spectrum()

    :param q_doas: queue.Queue   Queue where final processed dictionary is placed (should be a PyplisWorker.q_doas)
    """
    def __init__(self, routine=2, species=['SO2'], spec_specs=SpecSpecs(), spec_dir='C:\\', dark_dir=None,
                 q_doas=queue.Queue()):
        self.routine = routine          # Defines routine to be used, either (1) Polynomial or (2) Digital Filtering

        self.spec_specs = spec_specs    # Spectrometer specifications
        self.scan_proc = ScanProcess()  # Scanning object
        self.series = EmissionSeries()  # Emission rate timeseries object
        self.save_obj = SaveSpectra(self, self.scan_proc)   # Save object
        self.ss_str = self.spec_specs.file_ss.replace('{}', '')
        # ======================================================================================================================
        # Initial Definitions
        # ======================================================================================================================
        self.ppmm_conversion = 2.7e15   # convert absorption cross-section in cm2/molecule to ppm.m (MAY NEED TO CHANGE THIS TO A DICTIONARY AS THE CONVERSION MAY DIFFER FOR EACH SPECIES?)

        self.shift = 0                  # Shift of spectrum in number of pixels
        self.shift_used = 0             # Shift actually used in best fit plot
        self.shift_tol = 2              # Shift tolerance (will process data at multiple shifts defined by tolerance)
        self.stretch = 0                # Stretch of spectrum
        self.stretch_tol = 0            # As shift_tol but for stretch
        self.stretch_adjuster = 0.0001  # Factor to scale stretch (needed if different spectrometers have different pixel resolutions otherwise the stretch applied may be in too large or too small stages)
        self.stretch_resample = 100     # Number of points to resample the spectrum by during stretching
        self._start_stray_pix = None    # Pixel space stray light window definitions
        self._end_stray_pix = None
        self._start_stray_wave = 280    # Wavelength space stray light window definitions
        self._end_stray_wave = 290
        self._start_fit_pix = None  # Pixel space fitting window definitions
        self._end_fit_pix = None
        self._start_fit_wave = 305  # Wavelength space fitting window definitions
        self._end_fit_wave = 320
        self.fit_window = None      # Fitting window, determined by set_fit_window()
        self.fit_window_ref = None  # Placeholder for shifted fitting window for the reference spectrum
        self.wave_fit = True        # If True, wavelength parameters are used to define fitting window

        self.wavelengths = None         # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        self.wavelengths_cut = None     # Wavelengths in fit window
        self._dark_spec = None           # Dark spectrum
        self.dark_dict = {}             # Dictionary holding all dark spectra loaded in
        self._clear_spec_raw = None     # Clear (fraunhofer) spectrum - not dark corrected
        self._plume_spec_raw = None     # In-plume spectrum (main one which is used for calculation of SO2
        self.clear_spec_corr = None     # Clear (fraunhofer) spectrum - typically dark corrected and stray light corrected
        self.plume_spec_corr = None     # In-plume spectrum (main one which is used for calculation of SO2
        self.ref_spec = dict()          # Create empty dictionary for holding reference spectra
        self.ref_spec_interp = dict()   # Empty dictionary to hold reference spectra after sampling to spectrometer wavelengths
        self.ref_spec_conv = dict()     # Empty dictionary to hold reference spectra after convolving with ILS
        self.ref_spec_cut = dict()      # Ref spectrum cut to fit window
        self.ref_spec_ppmm = dict()   # Convolved ref spectrum scaled by ppmm_conversion factor
        self.ref_spec_filter = dict()   # Filtered reference spectrum
        self.ref_spec_fit = dict()      # Ref spectrum scaled by ppmm (for plotting)
        self.ref_spec_types = ['SO2', 'O3', 'ring'] # List of reference spectra types accepted/expected
        self.ref_spec_used  = species    # Reference spectra we actually want to use at this time (similar to ref_spec_types - perhaps one is obsolete (or should be!)
        self.abs_spec = None
        self.abs_spec_cut = None
        self.abs_spec_filt = None
        self.abs_spec_species = dict()  # Dictionary of absorbances isolated for individual species
        self.ILS_wavelengths = None     # Wavelengths for ILS
        self._ILS = None                 # Instrument line shape (will be a numpy array)
        self.processed_data = False     # Bool to define if object has processed DOAS yet - will become true once process_doas() is run
        self.auto_plume_params = False

        self.poly_order = 2  # Order of polynomial used to fit residual
        (self.filt_B, self.filt_A) = signal.butter(10, 0.065, btype='highpass')

        self.start_ca = -2000  # Starting column amount for iterations
        self.end_ca = 20000  # Ending column amount for iterations
        self.vals_ca = np.arange(self.start_ca, self.end_ca+1)  # Array of column amounts to be iterated over
        self.vals_ca_cut_idxs = np.arange(0, len(self.vals_ca), 100)
        self.vals_ca_cut = self.vals_ca[self.vals_ca_cut_idxs]
        self.mse_vals_cut = np.zeros(len(self.vals_ca_cut))
        self.mse_vals = np.zeros(len(self.vals_ca))  # Array to hold mse values

        self.std_err = None
        self.column_density = dict()

        self.filetypes = dict(defaultextension='.png', filetypes=[('PNG', '*.png')])

        # ----------------------------------------------------------------------------------
        # We need to make sure that all images are dark subtracted before final processing
        # Also make sure that we don't dark subtract more than once!
        self.ref_convolved = False  # Bool defining if reference spe has been convolved - speeds up DOAS processing
        self.new_spectra = True
        self.dark_corrected_clear = False
        self.dark_corrected_plume = False
        self.stray_corrected_clear = False    # Bool defining if all necessary spectra have been stray light corrected
        self.stray_corrected_plume = False    # Bool defining if all necessary spectra have been stray light corrected

        self.have_dark = False  # Used to define if a dark image is loaded.
        self.use_dark_dir = 0   # If true, a directory is used to hunt for dark spectra with correct shutter speed
        self.cal_dark_corr = False  # Tells us if the calibration image has been dark subtracted
        self.clear_dark_corr = False  # Tells us if the clear image has been dark subtracted
        self.plume_dark_corr = False  # Tells us if the plume image has been dark subtracted
        # ==============================================================================================================

        # Processing loop attributes
        self.process_thread = None      # Thread for running processing loop
        self.processing_in_thread = False   # Flags whether the object is processing in a thread or in the main thread - therefore deciding whether plots should be updated herein or through pyplisworker
        self.q_spec = queue.Queue()     # Queue where spectra files are placed, for processing herein
        self.q_doas = q_doas
        self.q_scan = queue.Queue()     # Queue where scan directories are placed for processing
        self.q_emission = queue.Queue() # Queue where emission rate and time are placed

        self._dark_dir = None
        self.dark_dir = dark_dir        # Directory where dark images are stored
        self.spec_dir = spec_dir        # Directory where plume spectra are stored
        self.date_dir = None            # Directory containing scan subdirectories. Should be date in format YYYY-MM-DD
        self.watch_dir = None           # Directory to watch, for automated program
        self.processed_dirs = []        # List of processed scan directories - to ensure we don't process again
        self.spec_dict = {}             # Dictionary containing all spectrum files from current spec_dir

        # Figures
        self.fig_spec = None            # plotting_gui.SpectraPlot object
        self.fig_doas = None            # plotting_gui.DOASPlot object
        self.fig_scan = None            # plotting_gui.CDPlot object
        self.fig_series = None          # plotting_gui.TimeSeriesPlot


    @property
    def start_stray_wave(self):
        return self._start_stray_wave

    @start_stray_wave.setter
    def start_stray_wave(self, value):
        self._start_stray_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._start_stray_pix = np.argmin(np.absolute(self.wavelengths - value))

    @property
    def end_stray_wave(self):
        return self._end_stray_wave

    @end_stray_wave.setter
    def end_stray_wave(self, value):
        self._end_stray_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._end_stray_pix = np.argmin(np.absolute(self.wavelengths - value))

    @property
    def start_fit_wave(self):
        return self._start_fit_wave

    @start_fit_wave.setter
    def start_fit_wave(self, value):
        self._start_fit_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._start_fit_pix = np.argmin(np.absolute(self.wavelengths - value))

    @property
    def end_fit_wave(self):
        return self._end_fit_wave

    @end_fit_wave.setter
    def end_fit_wave(self, value):
        self._end_fit_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._end_fit_pix = np.argmin(np.absolute(self.wavelengths - value))

    # --------------------------------------------
    # Set spectra attributes so that whenever they are updated we flag that they have not been dark or stray corrected
    # If we have a new dark spectrum too we need to assume it is to be used for correcting spectra, so all corrected
    # spectra, both dark and stray, become invalid
    @property
    def dark_spec(self):
        return self._dark_spec

    @dark_spec.setter
    def dark_spec(self, value):
        self._dark_spec = value

        # If we have a new dark image all dark and stray corrections become invalid
        self.dark_corrected_clear = False
        self.dark_corrected_plume = False
        self.stray_corrected_clear = False
        self.stray_corrected_plume = False

    @property
    def dark_dir(self):
        return self._dark_dir

    @dark_dir.setter
    def dark_dir(self, value):
        """If dark_dir is changed we need to reset the dark_dict which holds preloaded dark specs"""
        self.dark_dict = {}
        self._dark_dir = value

    @property
    def clear_spec_raw(self):
        return self._clear_spec_raw

    @clear_spec_raw.setter
    def clear_spec_raw(self, value):
        self._clear_spec_raw = value
        self.dark_corrected_clear = False
        self.stray_corrected_clear = False

    @property
    def plume_spec_raw(self):
        return self._plume_spec_raw

    @plume_spec_raw.setter
    def plume_spec_raw(self, value):
        self._plume_spec_raw = value
        self.dark_corrected_plume = False
        self.stray_corrected_plume = False

    @property
    def ILS(self):
        """ILS array"""
        return self._ILS

    @ILS.setter
    def ILS(self, value):
        self._ILS = value

        # If new ILS is generated, then must flag that ref spectrum is no longer convolved with up-to-date ILS
        self.ref_convolved = False
    # -------------------------------------------
    def set_doas_filename(self, filename):
        """Organise doas filename for saving purposes"""
        # Extract filename
        filename = os.path.split(filename)[-1]

        filename = filename.split('.')[0].replace('clear', 'doas')
        filename = filename.split('.')[0].replace('plume', 'doas')  # Filename may contain plume or clear
        doas_filename = filename + '_1.txt'
        self.doas_path = self.save_path + doas_filename

        idx = 2
        while os.path.exists(self.doas_path):
            self.doas_path = self.save_path + filename + '_{}.txt'.format(idx)
            idx += 1

    def get_spec_time(self, filename):
        """
        Gets time from filename and converts it to datetime object
        :param filename:
        :return spec_time:
        """
        # Make sure filename only contains file and not larger pathname
        filename = os.path.split(filename)[-1]

        # Extract time string from filename
        time_str = filename.split('_')[self.spec_specs.file_datestr_loc]

        # Turn time string into datetime object
        spec_time = datetime.datetime.strptime(time_str, self.spec_specs.file_datestr)

        return spec_time

    def dark_corr_spectra(self):
        """Subtract dark spectrum from spectra"""
        if not self.dark_corrected_clear:
            self.clear_spec_corr = self.clear_spec_raw - self.dark_spec
            self.clear_spec_corr[self.clear_spec_corr < 0] = 0
            self.dark_corrected_clear = True

        if not self.dark_corrected_plume:
            self.plume_spec_corr = self.plume_spec_raw - self.dark_spec
            self.plume_spec_corr[self.plume_spec_corr < 0] = 0
            self.dark_corrected_plume = True

    def stray_corr_spectra(self):
        """Correct spectra for stray light - spectra are assumed to be dark-corrected prior to running this function"""
        # Set the range of stray pixels
        stray_range = np.arange(self._start_stray_pix, self._end_stray_pix + 1)

        # Correct clear and plume spectra (assumed to already be dark subtracted)
        if not self.stray_corrected_clear:
            self.clear_spec_corr = self.clear_spec_corr - np.mean(self.clear_spec_corr[stray_range])
            self.clear_spec_corr[self.clear_spec_corr < 0] = 0
            self.stray_corrected_clear = True

        if not self.stray_corrected_plume:
            self.plume_spec_corr = self.plume_spec_corr - np.mean(self.plume_spec_corr[stray_range])
            self.plume_spec_corr[self.plume_spec_corr < 0] = 0
            self.stray_corrected_plume = True

    def load_ref_spec(self, pathname, species):
        """Load raw reference spectrum"""
        self.ref_spec[species] = np.loadtxt(pathname)

        # Remove un-needed data above 400 nm, which may slow down processing
        idxs = np.where(self.ref_spec[species][:, 0] > 400)
        self.ref_spec[species] = self.ref_spec[species][:np.min(idxs), :]

        # Assume we have loaded a new spectrum, so set this to False - ILS has not been convolved yet
        self.ref_convolved = False

    def get_ref_spectrum(self):
        """Load in reference spectrum"""
        self.wavelengths = None  # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        #
        # --------------------------------

    def conv_ref_spec(self):
        """
        Convolves reference spectrum with instument line shape (ILS)
        after first interpolating to spectrometer wavelengths
        """
        if self.wavelengths is None:
            print('No wavelength data to perform convolution')
            return

        # Need an odd sized array for convolution, so if even we omit the last pixel
        if self.ILS.size % 2 == 0:
            self.ILS = self.ILS[:-1]

        # Loop through all reference species we have loaded and resample their data to the spectrometers wavelengths
        for f in self.ref_spec_used:
            self.ref_spec_interp[f] = np.interp(self.wavelengths, self.ref_spec[f][:, 0], self.ref_spec[f][:, 1])
            self.ref_spec_conv[f] = convolve(self.ref_spec_interp[f], self.ILS)
            self.ref_spec_ppmm[f] = self.ref_spec_conv[f] * self.ppmm_conversion

        # Update bool as we have now performed this process
        self.ref_convolved = True

    def load_calibration_spectrum(self, pathname):
        """Load Calibation image for spectrometer"""
        pass

    def set_fit_windows(self):
        """Define fitting window for DOAS procedure
        If wavelength domain is used, first convert this to pixel space"""
        if self.wave_fit:
            if self.wavelengths is None:
                print('Error, first run get_ref_spectrum() to define wavelengths vector')
                return

        self.fit_window = np.arange(self._start_fit_pix, self._end_fit_pix)  # Fitting window (in Pixel space)
        self.fit_window_ref = self.fit_window - self.shift

    def stretch_spectrum(self, ref_key):
        """Stretch/squeeze reference spectrum to improve fit"""
        if self.stretch == 0:
            # If no stretch, extract reference spectrum using fit window and return
            return self.ref_spec_filter[ref_key][self.fit_window_ref]
        else:
            stretch_inc = (self.stretch * self.stretch_adjuster) / self.stretch_resample  # Stretch increment for resampled spectrum

            if self.stretch < 0:
                # Generate the fit window for extraction
                # Must be larger than fit window as a squeeze requires pulling in data outside of the fit window
                if self.fit_window_ref[-1] < (len(self.wavelengths) - 50):
                    extract_window = np.arange(self.fit_window_ref[0], self.fit_window_ref[-1] + 50)
                else:
                    extract_window = np.arange(self.fit_window_ref[0], len(self.wavelengths))
            else:
                extract_window = self.fit_window_ref

            wavelengths = self.wavelengths[extract_window]
            values = self.ref_spec_filter[ref_key][extract_window]

            # Generate new arrays with 'stretch_resample' more data points
            wavelengths_resampled = np.linspace(wavelengths[0], wavelengths[-1], len(wavelengths)*self.stretch_resample)
            values_resample = np.interp(wavelengths_resampled, wavelengths, values)

            num_pts = len(wavelengths_resampled)
            wavelengths_stretch = np.zeros(num_pts)

            # Stretch wavelengths
            for i in range(num_pts):
                wavelengths_stretch[i] = wavelengths_resampled[i] + (i * stretch_inc)

            values_stretch = np.interp(self.wavelengths[self.fit_window_ref], wavelengths_stretch, values_resample)

            return values_stretch

    def load_spec(self):
        """Load spectrum"""
        pass

    def load_dark(self):
        """Load drk images -> co-add to generate single dark image"""
        pass

    def reset_self(self):
        """
        Resets parameters associated with spectra, so that we can start from scratch without the potential for
        mistakenly using old data
        :return:
        """
        self.wavelengths = None
        self.plume_spec_raw = None
        self.clear_spec_raw = None

    def load_dir(self, prompt=True, plot=True):
        """Load spectrum directory
        :param: prompt  bool    If true, a dialogue box is opened to request directory load. Else, self.spec_dir is used
        :param: plot    bool    If true, the first spectra are plotted in the GUI"""

        if prompt:
            spec_dir = filedialog.askdirectory(title='Select spectrum sequence directory', initialdir=self.spec_dir)

            if len(spec_dir) > 0 and os.path.exists(spec_dir):
                self.spec_dir = spec_dir
            else:
                raise ValueError('Spectrum directory not recognised: {}'.format(spec_dir))
        else:
            if self.spec_dir is None:
                raise ValueError('Spectrum directory not recognised: {}'.format(self.spec_dir))

        # Update first_spec flag TODO possibly not used in DOASWorker, check
        self.first_spec = True

        # Reset self
        self.reset_self()

        # Get list of all files in directory
        self.spec_dict = self.get_spec_list()

        if len(self.spec_dict['plume']) == 0:
            print('No plume files present in directory')
            return

        # Set current spectra to first in lists
        if len(self.spec_dict['clear']) > 0:
            self.wavelengths, self.clear_spec_raw = load_spectrum(self.spec_dir + self.spec_dict['clear'][0])
            self.clear_path = self.spec_dir + self.spec_dict['clear'][0]
        else:
            # If no clear spectra are in the directory, we use teh first plume image as the clear sky
            self.wavelengths, self.clear_spec_raw = load_spectrum(self.spec_dir + self.spec_dict['plume'][0])
            self.clear_path = self.spec_dir + self.spec_dict['plume'][0]
        if len(self.spec_dict['plume']) > 0:
            self.wavelengths, self.plume_spec_raw = load_spectrum(self.spec_dir + self.spec_dict['plume'][0])
            self.plume_path = self.spec_dir + self.spec_dict['plume'][0]
        if len(self.spec_dict['dark']) > 0:
            ss_id = self.spec_specs.file_ss.replace('{}', '')
            ss = self.spec_dict['plume'][0].split('_')[self.spec_specs.file_ss_loc].replace(ss_id, '')

            # If filename is not as expected (we don't have the ss flag) we set ss to None. Any dark spectrum from dark_dir
            # will then be accepted as teh dark spec
            if ss_id not in ss:
                ss = None
            self.dark_spec = self.find_dark_spectrum(self.dark_dir, ss)

        # Update plots if requested
        if plot:
            self.fig_spec.update_clear()
            self.fig_spec.update_dark()
            self.fig_spec.update_plume()

    def get_spec_list(self):
        """
        Gets list of spectra files from current spec_dir
        Assumes all .npy files in the directory are spectra...
        :returns: sd    dict    Dictionary containing all filenames
        """
        # Setup empty dictionary sd
        sd = {}

        # Get all files into associated list/dictionary entry (search for npy or txt files)
        sd['all'] = [f for f in os.listdir(self.spec_dir) if self.spec_specs.file_ext in f or '.txt' in f]
        sd['all'].sort()
        sd['plume'] = [f for f in sd['all']
                       if self.spec_specs.file_spec_type['meas'].lower() in f.lower()]
        sd['plume'].sort()
        sd['clear'] = [f for f in sd['all']
                       if self.spec_specs.file_spec_type['clear'].lower() in f.lower()]
        sd['clear'].sort()
        sd['dark'] = [f for f in sd['all']
                      if self.spec_specs.file_spec_type['dark'].lower() in f.lower()]
        sd['dark'].sort()
        return sd

    def get_ss_from_filename(self, filename):
        """
        Extracts shutter speeds from filename
        :param filename:
        :return ss: str Shutter speed
        """
        ss_full_str = filename.split('_')[self.spec_specs.file_ss_loc]
        if self.ss_str not in ss_full_str:
            ss = None
        else:
            ss = int(ss_full_str.replace(self.ss_str, ''))
        return ss

    def find_dark_spectrum(self, spec_dir, ss):
        """
        Searches for suitable dark spectrum in designated directory by finding one with the same shutter speed as
        passed to function.
        :return: dark_spec
        """
        # If ss is None, we get the first dark spectrum in the directory
        if ss is None:
            # List all dark images in directory
            dark_list = [f for f in os.listdir(spec_dir)
                         if self.spec_specs.file_spec_type['dark'].lower() in f.lower()
                         and (self.spec_specs.file_ext in f or '.txt' in f)]
            if len(dark_list) > 0:
                dark_spec = load_spectrum(spec_dir + dark_list[0])[1]
            else:
                dark_spec = None
            return dark_spec

        # Ensure ss is a string
        ss = str(ss)

        # Fast dictionary look up for preloaded dark spectra
        if ss in self.dark_dict.keys():
            dark_spec = self.dark_dict[ss]
            return dark_spec

        # List all dark images in directory
        dark_list = [f for f in os.listdir(spec_dir)
                     if self.spec_specs.file_spec_type['dark'].lower() in f.lower() and self.spec_specs.file_ext in f]

        # Extract ss from each image and round to 2 significant figures
        ss_list = [int(f.split('_')[self.spec_specs.file_ss_loc].replace(self.ss_str, '')) for f in dark_list]

        ss_idx = [i for i, x in enumerate(ss_list) if x == int(ss)]
        ss_spectra = [dark_list[i] for i in ss_idx]

        if len(ss_spectra) < 1:
            return None

        # If we have images, we loop through them to create a coadded image
        dark_full = np.zeros([self.spec_specs.pix_num, len(ss_spectra)])
        for i, ss_spectrum in enumerate(ss_spectra):
            # Load image. Coadd.
            wavelengths, dark_full[:, i] = load_spectrum(spec_dir + ss_spectrum)

        # Coadd images to creat single image
        dark_spec = np.mean(dark_full, axis=1)

        # Update lookup dictionary for fast retrieval of dark image later
        self.dark_dict[ss] = dark_spec

        return dark_spec

    def save_dark(self, filename):
        """Save dark spectrum"""
        if self.wavelengths is None or self.dark_spec is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.dark_spec):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.dark_spec]),
                   header='Dark spectrum\nWavelength [nm]\tIntensity [DN]')

    def save_clear_raw(self, filename):
        """Save clear spectrum"""
        if self.wavelengths is None or self.clear_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.clear_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.clear_spec_raw]),
                   header='Raw clear (Fraunhofer) spectrum\n'
                          '-Not dark-corrected\nWavelength [nm]\tIntensity [DN]')

    def save_plume_raw(self, filename):
        if self.wavelengths is None or self.plume_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.plume_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.plume_spec_raw]),
                   header='Raw in-plume spectrum\n'
                          '-Not dark-corrected\nWavelength [nm]\tIntensity [DN]')

    def load_plume_params(self):
        """Loads plume speed and plume distance parameters from file (for automated processing)"""
        # Initialise parameters in case they don't exist in the file
        plume_speed = None
        plume_distance = None

        # Plume file must be located in the scanning directory
        pathname = self.spec_dir + self.spec_specs.plume_params_file

        # Loop through file and extract plume speed and distance information
        with open(pathname, 'r') as f:
            for line in f:
                if self.spec_specs.plume_speed_id in line:
                    plume_speed = float(line.split(self.spec_specs.plume_speed_id)[1].split('\n')[0])
                elif self.spec_specs.plume_dist_id in line:
                    plume_distance = float(line.split(self.spec_specs.plume_dist_id)[1].split('\n')[0])

        return plume_speed, plume_distance

    @staticmethod
    def doas_fit(ref_spec, *fit_params):
        """
        DOAS fit function

        Parameters
        -----------
        ref_spec : array-like object
            Contains all (k) reference spectra to be fitted in a (k,N) array.
            N is the size of the spectrum grid
        """
        # Unpack fit scalars into column vector
        scalars = np.transpose(np.array([list(fit_params)]))

        # Sum reference spectra to create total absorbance for comparison with measured absorbance spectrum
        total_absorbance = np.sum(ref_spec * scalars, axis=0)

        return total_absorbance

    def poly_doas(self):
        """
        Performs main processing in polynomial fitting DOAS retrieval
        NOT COMPLETE, DO NOT USE!!!
        """
        self.wavelengths_cut = self.wavelengths[self.fit_window]  # Extract wavelengths (used in plotting)

        with np.errstate(divide='ignore'):
            self.abs_spec = np.log(np.divide(self.clear_spec_corr, self.plume_spec_corr))  # Calculate absorbance

        self.ref_spec_cut['SO2'] = self.ref_spec_ppmm[self.ref_spec_types['SO2']][self.fit_window_ref]
        self.abs_spec_cut = self.abs_spec[self.fit_window]

        idx = 0
        for i in self.vals_ca:
            ref_spec_fit = self.ref_spec_scaled['SO2'] * i  # Our iterative guess at the SO2 column density
            residual = self.abs_spec_cut - ref_spec_fit  # Calculate resultant residual from spectrum fitting
            poly_fit = np.polyfit(self.fit_window, residual, self.poly_order)  # Fit polynomial to residual
            poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window

            self.mse_vals[idx] = np.mean(np.power(residual - poly_vals, 2))  # Calculate MSE of fit

            idx += 1

        self.min_idx = np.argmin(self.mse_vals)
        self.column_density['SO2'] = self.vals_ca[self.min_idx]

    def fltr_doas(self):
        """
        Performs main retrieval in digital filtering DOAS retrieval
        """
        self.wavelengths_cut = self.wavelengths[self.fit_window]    # Extract wavelengths (used in plotting)

        # Calculate absorbance
        with np.errstate(divide='ignore', invalid='ignore'):
            self.abs_spec = np.log(np.divide(self.clear_spec_corr, self.plume_spec_corr))

        # Remove nans and infs (filter function doesn't work if they are included)
        self.abs_spec[np.isinf(self.abs_spec)] = np.nan
        self.abs_spec[np.isnan(self.abs_spec)] = 0

        self.abs_spec_filt = signal.lfilter(self.filt_B, self.filt_A, self.abs_spec)    # Filter absorbance spectrum
        self.abs_spec_cut = self.abs_spec_filt[self.fit_window]                         # Extract fit window

        # Loop through reference spectra and prepare them for processing
        for spec in self.ref_spec_used:

            # Filter reference spectrum and extract fit window
            self.ref_spec_filter[spec] = signal.lfilter(self.filt_B, self.filt_A, self.ref_spec_ppmm[spec])

            # Stretch spectrum
            self.ref_spec_cut[spec] = self.stretch_spectrum(spec)

        # self.ref_spec_cut['SO2'] = self.ref_spec_filter['SO2'][self.fit_window_ref]

        # ------------------------------------------------------------------------------------------
        # Scipy.optimize.curve_fit least squares fitting
        # ------------------------------------------------------------------------------------------
        # Pack all requested reference spectra into an array for curve fitting
        ref_spectra_packed = np.empty((len(self.ref_spec_used), len(self.abs_spec_cut)))
        i = 0
        for spec in self.ref_spec_used:
            ref_spectra_packed[i, :] = self.ref_spec_cut[spec]
            i += 1

        # Run fit
        column_densities, pcov = curve_fit(self.doas_fit, ref_spectra_packed, self.abs_spec_cut,
                                           p0=np.ones(ref_spectra_packed.shape[0]))
        self.std_err = round(np.sqrt(np.diag(pcov))[0], 1)

        # Loop through species to unpack results
        i = 0
        self.ref_spec_fit['Total'] = np.zeros(len(self.abs_spec_cut))
        for spec in self.ref_spec_used:
            # Unpack column densities to dictionary
            self.column_density[spec] = int(round(column_densities[i]))

            # Generate scaled reference spectrum
            self.ref_spec_fit[spec] = self.ref_spec_cut[spec] * self.column_density[spec]

            # Add scaled reference spectrum to total ref spec
            self.ref_spec_fit['Total'] += self.ref_spec_fit[spec]

            i += 1

    def gen_abs_specs(self):
        """
        Generate absorbance spectra for individual species and match to measured absorbance
        :return:
        """
        # Iterate through reference spectra
        for spec in self.ref_spec_used:

            # Retrieve all ref spectra other than that held in spec
            other_species = [f for f in self.ref_spec_used if f is not spec]

            # Subtract contributions of these other spectra to absorption
            self.abs_spec_species[spec] = self.abs_spec_cut
            for i in other_species:
                self.abs_spec_species[spec] = self.abs_spec_species[spec] - self.ref_spec_fit[i]

        # Setup total absorbance spectrum as 'Total' key
        self.abs_spec_species['Total'] = self.abs_spec_cut

        # Calculate final residual by subtracting all absorbing species
        self.abs_spec_species['residual'] = self.abs_spec_cut
        for spec in self.ref_spec_used:
            self.abs_spec_species['residual'] = self.abs_spec_species['residual'] - self.ref_spec_fit[spec]


    def poly_plot_gen(self):
        """Generate arrays to be plotted -> residual, fitted spectrum"""
        self.ref_spec_fit['SO2'] = self.ref_spec_cut['SO2'] * self.column_density['SO2']
        self.residual = self.abs_spec_cut - self.ref_spec_fit['SO2']
        poly_fit = np.polyfit(self.fit_window, self.residual, self.poly_order)  # Fit polynomial to residual
        self.poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window
        self.best_fit = self.ref_spec_fit['SO2'] + self.poly_vals  # Generate best fit absorbance spectrum

        # MAKE PLOT
        plt.figure()
        abs_plt, = plt.plot(self.abs_spec_cut, label='Absorbance spectrum')
        ref_plt, = plt.plot(self.ref_spec_fit, label='Reference spectrum * CA')
        res_plt, = plt.plot(self.residual, label='Residual')
        poly_plt, = plt.plot(self.poly_vals, label='Polynomial fit')
        best_plt, = plt.plot(self.best_fit, label='Best fit')
        plt.xlabel('Pixel')
        plt.ylabel('Absorbance')
        plt.legend(handles=[abs_plt, ref_plt, res_plt, poly_plt, best_plt])
        plt.show()

    def process_doas(self):
        """Handles the order of DOAS processing"""
        # Check we have all of the correct spectra to perform processing
        if self.clear_spec_raw is None or self.plume_spec_raw is None or self.wavelengths is None:
            raise SpectraError('Require clear and plume spectra for DOAS processing')

        if self.ref_spec_types[0] not in self.ref_spec.keys():
            raise SpectraError('No SO2 reference spectrum present for processing')

        # If the stray region hasn't been converted to pixel space we just set it to itself to implement pixel mapping
        if self._start_stray_pix is None or self._end_stray_pix is None:
            self.start_stray_wave = self.start_stray_wave
            self.end_stray_wave = self.end_stray_wave

        # Same for fit window
        if self._start_fit_pix is None or self._end_fit_pix is None:
            self.start_fit_wave = self.start_fit_wave
            self.end_fit_wave = self.end_fit_wave

        if not self.dark_corrected_clear or not self.dark_corrected_plume:
            if self.dark_spec is None:
                print('Warning! No dark spectrum present, processing without dark subtraction')

                # Set raw spectra to the corrected spectra, ignoring that they have not been dark corrected
                self.clear_spec_corr = self.clear_spec_raw
                self.plume_spec_corr = self.plume_spec_raw
            else:
                self.dark_corr_spectra()

        # Correct spectra for stray light
        if not self.stray_corrected_clear or not self.stray_corrected_plume:
            self.stray_corr_spectra()

        # Convolve reference spectrum with the instrument lineshape
        # TODO COnfirm this still works - prior to the shift_tol incorporation these lines came AFTER self.set_fit_windows()
        # TODO I think it won't make a difference as they aren't directly related. But I should confirm this
        if not self.ref_convolved:
            self.conv_ref_spec()

        # Process at range of shifts defined by shift tolerance
        fit_results = []
        for i, shift in enumerate(range(self.shift - self.shift_tol, self.shift + self.shift_tol + 1)):
            # Set shift to new value
            self.shift = shift

            # Set fitting windows for acquired and reference spectra
            self.set_fit_windows()

            # Run processing
            self.fltr_doas()

            # Create dictionary of all relevant fit results (need to make copies, otherwise dictionaries change when changed elsewhere)
            fit_results.append({'std_err': self.std_err,
                                'column_density': copy.deepcopy(self.column_density),
                                'ref_spec_fit': copy.deepcopy(self.ref_spec_fit),
                                'shift_used': shift
                                })

        # Find best fit by looping through all results, extracting std_err and finding minimum value
        best_fit = np.argmin(np.array([x['std_err'] for x in fit_results]))
        best_fit_results = fit_results[best_fit]

        # Unpack the fit results dictionary back into the object attributes the keys represent
        for key in best_fit_results:
            setattr(self, key, best_fit_results[key])

        # Reset shift to starting value (otherwise we could possibly get a runaway shift during multiple processing?)
        self.shift -= self.shift_tol

        # Generate absorbance spectra for individual species
        self.gen_abs_specs()

        # Set flag defining that data has been fully processed
        self.processed_data = True

    def process_scan(self, scan_dir=None, plume_speed=None, plume_distance=None, plot=True):
        """
        Processes a scan to calculate emission rates
        :param scan_dir: str            Directory containing scan data
        :param plume_speed: float       Plume speed [m/s]
        :param plume_distance: float    Distance to plume from scanner
        :param plot: bool               Defines whether plots should be updated
        :return:
        """
        # Update parameters if requested
        if scan_dir is not None:
            self.spec_dir = scan_dir

        # If auto plume params are set then we attempt to load them
        if self.auto_plume_params:
            try:
                self.scan_proc.plume_speed, self.scan_proc.plume_distance = self.load_plume_params()
                # If the values after attempting to load from file are None, then we must take the values from the
                # GUI. Otherwise we use the loaded values.
                if self.scan_proc.plume_speed is None:
                    self.scan_proc.plume_speed = self.fig_scan.plume_speed
                else:
                    self.fig_scan.plume_speed = self.scan_proc.plume_speed
                if self.scan_proc.plume_distance is None:
                    self.scan_proc.plume_distance = self.fig_scan.plume_dist
                else:
                    self.fig_scan.plume_dist = self.scan_proc.plume_distance
            # If an exception is raised when trying to load from file, we take the values from the GUI
            except BaseException:
                # Update scan parameters
                self.scan_proc.plume_distance = self.fig_scan.plume_dist
                self.scan_proc.plume_speed = self.fig_scan.plume_speed
        else:
            # Update scan parameters
            self.scan_proc.plume_distance = self.fig_scan.plume_dist
            self.scan_proc.plume_speed = self.fig_scan.plume_speed

            if plume_speed is not None:
                self.scan_proc.plume_speed = plume_speed

            if plume_distance is not None:
                self.scan_proc.plume_distance = plume_distance

        # Dark spectrum should be held in same directory
        if not self.use_dark_dir:
            self.dark_dir = self.spec_dir

        # Create directory to save processed data
        self.save_path = self.spec_dir + 'Processing_1/'
        i = 2
        while os.path.exists(self.save_path):
            self.save_path = self.spec_dir + 'Processing_{}/'.format(i)
            i += 1
        os.mkdir(self.save_path)

        # Initialise values
        self.scan_proc.clear_data()
        scan_angle = 0

        # Load directory - this loads in clear spectrum, else uses first plume spectrum as clear. Generates file list
        # for each spectrum type. Sets up dark image from scan
        self.load_dir(prompt=False, plot=plot)

        # Update dark and spectrum, which will remain the same through the scan processing
        self.fig_spec.update_dark()
        self.fig_spec.update_clear()
        self.fig_scan.clear_plot()

        # Loop through plume files and process them
        for plume_file in self.spec_dict['plume']:
            # Setup full path
            pathname = self.spec_dir + plume_file

            # Extract datetime object of spectrum from filename
            spec_time = self.get_spec_time(plume_file)

            # ----------------------------------------------------------------------
            # USING DARK DIRECTORY UPDATES
            # Extract shutter speed
            if self.use_dark_dir:
                filename = os.path.split(pathname)[-1]
                ss = self.get_ss_from_filename(filename)
                # Find dark spectrum with same shutter speed
                self.dark_spec = self.find_dark_spectrum(self.dark_dir, ss)
            # ---------------------------------------------------------------------

            # Load spectrum and update figure
            self.wavelengths, self.plume_spec_raw = load_spectrum(pathname)
            # ===========================================
            # EDIT FOR WAVELENGTHS BEING IN MILLIONS FOR IGP
            self.wavelengths = self.wavelengths/10000
            # ===========================================
            self.fig_spec.update_plume()

            # Process spectrum and update plot
            self.process_doas()
            self.fig_doas.update_plot()

            # Gather all relevant information and spectra
            processed_dict = {'processed': True,  # Flag whether this is a complete, processed dictionary
                              'time': spec_time,
                              'scan_angle': scan_angle, # Scan angle
                              'filename': pathname,  # Filename of processed spectrum
                              'dark': self.dark_spec,  # Dark spectrum used (raw)
                              'clear': self.clear_spec_raw,  # Clear spectrum used (raw)
                              'plume': self.plume_spec_raw,  # Plume spectrum used (raw)
                              'shift': self.shift_used,     # Shift used in fit - may depend on shift tolerance
                              'abs': self.abs_spec_cut,  # Absorption spectrum
                              'ref': self.abs_spec_species,  # Reference spectra (scaled)
                              'std_err': self.std_err,  # Standard error of fit
                              'column_density': self.column_density}  # Column density

            # Update scan processor
            self.scan_proc.add_data(processed_dict['scan_angle'], processed_dict['column_density']['SO2'])
            self.fig_scan.update_plot()

            # Put results in doas q to be plotted
            self.q_doas.put(processed_dict)

            # Save processed data
            self.set_doas_filename(pathname)
            self.save_obj.save_processed_spec('Scan directory dark', self.clear_path, pathname, self.doas_path)

            # Step the scan angle
            scan_angle += self.scan_proc.scan_incr

        # Calculate emission rate
        self.scan_proc.calc_emission_rate()
        self.fig_scan.update_emission_rate()

        # Save scan data
        self.save_obj.save_scan(self.save_path)

        # Put this emission rate in q with scan dir
        self.q_emission.put([processed_dict['time'], self.scan_proc.SO2_flux])
        self.series.add_data(processed_dict['time'], self.scan_proc.SO2_flux)
        self.fig_series.update_plot()
        self.fig_series.update_emission_rate()

    def start_continuous_processing(self):
        """
        Starts continuous processing of scan directories
        :return:
        """
        self.processing_in_thread = True
        self.process_thread = threading.Thread(target=self._process_scan_loop, args=())
        self.process_thread.daemon = True
        self.process_thread.start()

        # Reset time series
        self.series.clear_data()

        # Setup watcher thread to watch directory and then add new directories to q_scan
        # Setup directory watcher. recursive=True, to watch subdirectories
        if self.watch_dir is None:
            messagebox.showerror('No directory selected',
                                 'Please select a directory to watch before attempting to start the directory watcher')
            return
        self.scan_dir_watcher = create_dir_watcher(self.watch_dir, True, self.directory_watch_handler)
        self.scan_dir_watcher.start()

    def stop_continuous_processing(self):
        """Stops directory watcher and continuous processing"""
        self.scan_dir_watcher.stop()
        self.q_scan.put('exit')

    def directory_watch_handler(self, pathname, t):
        """Controls the watching of a directory"""
        print('Got file: {}'.format(pathname))

        # Separate the filename and pathname
        pathname, filename = os.path.split(pathname)

        if 'Processing_' in pathname:
            return
        if filename == self.series.filename:
            return

        # Only process this directory once the scan_complete file is present
        if filename == self.spec_specs.scan_complete:
            # Extract date directory
            date_path = os.path.split(pathname)[0]
            date_dir = os.path.split(date_path)[-1]

            # If this is data from a new day we must reset data
            if date_dir != self.date_dir:
                # Save time series for this day
                self.series.save_series(date_path, self.scan_proc.plume_distance, self.scan_proc.plume_speed)

                # Clear time series data
                self.series.clear_data()

            # Extract directory from full path
            scan_dir = pathname + '/'

            # Put directory in queue to be processed
            self.q_scan.put(scan_dir)

    def _process_scan_loop(self):
        """
        Main process loop for doas processing of scan directories
        :return:
        """
        while True:
            # Get scan directory
            scan_dir = self.q_scan.get(block=True)

            # Can end function by adding exit to queue
            if scan_dir == 'exit':
                return

            # Process scan
            self.process_scan(scan_dir=scan_dir, plot=False)

    def start_processing_threadless(self):
        """
        Process spectra already in a directory, without entering a thread - this means that the _process_loop
        function can update the plots wihtout going through the main thread in PyplisWorker
        """
        # Flag that we are running processing outside of thread
        self.processing_in_thread = False

        # Get all files
        spec_files = [f for f in os.listdir(self.spec_dir) if '.npy' in f]

        # Sort files alphabetically (which will sort them by time due to file format)
        spec_files.sort()

        # Extract clear spectra if they exist. If not, the first file is assumed to be the clear spectrum
        clear_spec = [f for f in spec_files if self.spec_specs.file_spec_type['clear'].lower() + '.npy' in f.lower()]
        plume_spec = [f for f in spec_files if self.spec_specs.file_spec_type['meas'].lower() + '.npy' in f.lower()]

        # Loop through all files and add them to queue
        for file in clear_spec:
            self.q_spec.put(self.spec_dir + file)
        for file in plume_spec:
            self.q_spec.put(self.spec_dir + file)

        # Add the exit flag at the end, to ensure that the process_loop doesn't get stuck waiting on the queue forever
        self.q_spec.put('exit')

        # Begin processing
        self._process_loop()


    def start_processing_thread(self):
        """Public access thread starter for _processing"""
        self.processing_in_thread = True
        self.process_thread = threading.Thread(target=self._process_loop, args=())
        self.process_thread.daemon = True
        self.process_thread.start()

    def _process_loop(self):
        """
        Main process loop for doas
        :return:
        """
        first_spec = True       # First spectrum is used as clear spectrum

        while True:
            # Blocking wait for new file
            pathname = self.q_spec.get(block=True)

            # Close thread if requested with 'exit' command
            if pathname == 'exit':
                break

            # Extract filename and create datetime object of spectrum time
            filename = os.path.split(pathname)[-1]
            spec_time = self.get_spec_time(filename)

            # Extract shutter speed
            ss_full_str = filename.split('_')[self.spec_specs.file_ss_loc]
            if self.ss_str not in ss_full_str:
                ss = None
            else:
                ss = int(ss_full_str.replace(self.ss_str, ''))

            # Find dark spectrum with same shutter speed
            self.dark_spec = self.find_dark_spectrum(self.dark_dir, ss)

            # Load spectrum
            self.wavelengths, spectrum = load_spectrum(pathname)

            # Make first spectrum clear spectrum
            if first_spec:
                self.clear_spec_raw = spectrum

                processed_dict = {'processed': False,
                                  'time': spec_time,
                                  'filename': pathname,
                                  'dark': self.dark_spec,
                                  'clear': self.clear_spec_raw}



            # Process plume spectrum
            else:
                self.plume_spec_raw = spectrum
                self.process_doas()

                # Gather all relevant information and spectra and pass it to PyplisWorker
                processed_dict = {'processed': True,             # Flag whether this is a complete, processed dictionary
                                  'time': spec_time,
                                  'filename': pathname,             # Filename of processed spectrum
                                  'dark': self.dark_spec,           # Dark spectrum used (raw)
                                  'clear': self.clear_spec_raw,     # Clear spectrum used (raw)
                                  'plume': self.plume_spec_raw,     # Plume spectrum used (raw)
                                  'abs': self.abs_spec_cut,         # Absorption spectrum
                                  'ref': self.abs_spec_species,     # Reference spectra (scaled)
                                  'std_err': self.std_err,          # Standard error of fit
                                  'column_density': self.column_density}    # Column density

            # Pass data dictionary to PyplisWorker queue
            if self.processing_in_thread:
                self.q_doas.put(processed_dict)

            # Or if we are not in a processing thread we update the plots directly here
            else:
                self.fig_spec.update_dark()
                if first_spec:
                    self.fig_spec.update_clear()
                else:
                    self.fig_spec.update_plume()

                    # Update doas plot
                    self.fig_doas.update_plot()

            if first_spec:
                # Now that we have processed first_spec, set flag to False
                first_spec = False


class DOASWorkerOld:
    """
    Old version of DOASWorker - used for SpecScan, but should now be becoming obsolete as I transfer improved version
    from pycam

    Class to control DOAS processing
    General order of play for processing:
    Initiate class,
    get_ref_spectrum()
    set_fit_window()
    shift_spectrum()
    """
    def __init__(self, routine=2, species=['SO2']):
        self.routine = routine  # Defines routine to be used, either (1) Polynomial or (2) Digital Filtering

        self.spec_specs = SpecSpecs()
        # ======================================================================================================================
        # Initial Definitions
        # ======================================================================================================================
        self.ppmm_conversion = 2.7e15   # convert absorption cross-section in cm2/molecule to ppm.m (MAY NEED TO CHANGE THIS TO A DICTIONARY AS THE CONVERSION MAY DIFFER FOR EACH SPECIES?)

        self.shift = 0                  # Shift of spectrum in number of pixels
        self.shift_tol = 2
        self.stretch = 0                # Stretch of spectrum
        self.stretch_tol = 0
        self.stretch_adjuster = 0.0001  # Factor to scale stretch (needed if different spectrometers have different pixel resolutions otherwise the stretch applied may be in too large or too small stages)
        self.stretch_resample = 100     # Number of points to resample the spectrum by during stretching
        self._start_stray_pix = None    # Pixel space stray light window definitions
        self._end_stray_pix = None
        self._start_stray_wave = 280    # Wavelength space stray light window definitions
        self._end_stray_wave = 290
        self._start_fit_pix = None  # Pixel space fitting window definitions
        self._end_fit_pix = None
        self._start_fit_wave = 305  # Wavelength space fitting window definitions
        self._end_fit_wave = 320
        self.fit_window = None      # Fitting window, determined by set_fit_window()
        self.fit_window_ref = None  # Placeholder for shifted fitting window for the reference spectrum
        self.wave_fit = True        # If True, wavelength parameters are used to define fitting window

        self.wavelengths = None         # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        self.wavelengths_cut = None     # Wavelengths in fit window
        self._dark_spec = None           # Dark spectrum
        self._clear_spec_raw = None     # Clear (fraunhofer) spectrum - not dark corrected
        self._plume_spec_raw = None     # In-plume spectrum (main one which is used for calculation of SO2
        self.clear_spec_corr = None     # Clear (fraunhofer) spectrum - typically dark corrected and stray light corrected
        self.plume_spec_corr = None     # In-plume spectrum (main one which is used for calculation of SO2
        self.ref_spec = dict()          # Create empty dictionary for holding reference spectra
        self.ref_spec_interp = dict()   # Empty dictionary to hold reference spectra after sampling to spectrometer wavelengths
        self.ref_spec_conv = dict()     # Empty dictionary to hold reference spectra after convolving with ILS
        self.ref_spec_cut = dict()      # Ref spectrum cut to fit window
        self.ref_spec_ppmm = dict()   # Convolved ref spectrum scaled by ppmm_conversion factor
        self.ref_spec_filter = dict()   # Filtered reference spectrum
        self.ref_spec_fit = dict()      # Ref spectrum scaled by ppmm (for plotting)
        self.ref_spec_types = ['SO2', 'O3', 'ring'] # List of reference spectra types accepted/expected
        self.ref_spec_used  = species    # Reference spectra we actually want to use at this time (similar to ref_spec_types - perhaps one is obsolete (or should be!)
        self.abs_spec = None
        self.abs_spec_cut = None
        self.abs_spec_filt = None
        self.abs_spec_species = dict()  # Dictionary of absorbances isolated for individual species
        self.ILS = None                 # Instrument line shape (will be a numpy array)
        self.processed_data = False     # Bool to define if object has processed DOAS yet - will become true once process_doas() is run

        self.poly_order = 2  # Order of polynomial used to fit residual
        (self.filt_B, self.filt_A) = signal.butter(10, 0.065, btype='highpass')

        self.start_ca = -2000  # Starting column amount for iterations
        self.end_ca = 20000  # Ending column amount for iterations
        self.vals_ca = np.arange(self.start_ca, self.end_ca+1)  # Array of column amounts to be iterated over
        self.vals_ca_cut_idxs = np.arange(0, len(self.vals_ca), 100)
        self.vals_ca_cut = self.vals_ca[self.vals_ca_cut_idxs]
        self.mse_vals_cut = np.zeros(len(self.vals_ca_cut))
        self.mse_vals = np.zeros(len(self.vals_ca))  # Array to hold mse values

        self.std_err = None
        self.column_density = dict()

        self.filetypes = dict(defaultextension='.png', filetypes=[('PNG', '*.png')])

        # ----------------------------------------------------------------------------------
        # We need to make sure that all images are dark subtracted before final processing
        # Also make sure that we don't dark subtract more than once!
        self.ref_convolved = False  # Bool defining if reference spe has been convolved - speeds up DOAS processing
        self.new_spectra = True
        self.dark_corrected_clear = False
        self.dark_corrected_plume = False
        self.stray_corrected_clear = False    # Bool defining if all necessary spectra have been stray light corrected
        self.stray_corrected_plume = False    # Bool defining if all necessary spectra have been stray light corrected

        self.have_dark = False  # Used to define if a dark image is loaded.
        self.cal_dark_corr = False  # Tells us if the calibration image has been dark subtracted
        self.clear_dark_corr = False  # Tells us if the clear image has been dark subtracted
        self.plume_dark_corr = False  # Tells us if the plume image has been dark subtracted
        # ==============================================================================================================

    @property
    def start_stray_wave(self):
        return self._start_stray_wave

    @start_stray_wave.setter
    def start_stray_wave(self, value):
        self._start_stray_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._start_stray_pix = np.argmin(np.absolute(self.wavelengths - value))

    @property
    def end_stray_wave(self):
        return self._end_stray_wave

    @end_stray_wave.setter
    def end_stray_wave(self, value):
        self._end_stray_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._end_stray_pix = np.argmin(np.absolute(self.wavelengths - value))


    @property
    def start_fit_wave(self):
        return self._start_fit_wave

    @start_fit_wave.setter
    def start_fit_wave(self, value):
        self._start_fit_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._start_fit_pix = np.argmin(np.absolute(self.wavelengths - value))

    @property
    def end_fit_wave(self):
        return self._end_fit_wave

    @end_fit_wave.setter
    def end_fit_wave(self, value):
        self._end_fit_wave = value

        # Set pixel value too, if wavelengths attribute is present
        if self.wavelengths is not None:
            self._end_fit_pix = np.argmin(np.absolute(self.wavelengths - value))

    # --------------------------------------------
    # Set spectra attributes so that whenever they are updated we flag that they have not been dark or stray corrected
    # If we have a new dark spectrum too we need to assume it is to be used for correcting spectra, so all corrected
    # spectra, both dark and stray, become invalid
    @property
    def dark_spec(self):
        return self._dark_spec

    @dark_spec.setter
    def dark_spec(self, value):
        self._dark_spec = value

        # If we have a new dark image all dark and stray corrections become invalid
        self.dark_corrected_clear = False
        self.dark_corrected_plume = False
        self.stray_corrected_clear = False
        self.stray_corrected_plume = False

    @property
    def clear_spec_raw(self):
        return self._clear_spec_raw

    @clear_spec_raw.setter
    def clear_spec_raw(self, value):
        self._clear_spec_raw = value
        self.dark_corrected_clear = False
        self.stray_corrected_clear = False

    @property
    def plume_spec_raw(self):
        return self._plume_spec_raw

    @plume_spec_raw.setter
    def plume_spec_raw(self, value):
        self._plume_spec_raw = value
        self.dark_corrected_plume = False
        self.stray_corrected_plume = False
    # -------------------------------------------

    def dark_corr_spectra(self):
        """Subtract dark spectrum from spectra"""
        if not self.dark_corrected_clear:
            self.clear_spec_corr = self.clear_spec_raw - self.dark_spec
            self.clear_spec_corr[self.clear_spec_corr < 0] = 0
            self.dark_corrected_clear = True

        if not self.dark_corrected_plume:
            self.plume_spec_corr = self.plume_spec_raw - self.dark_spec
            self.plume_spec_corr[self.plume_spec_corr < 0] = 0
            self.dark_corrected_plume = True

    def stray_corr_spectra(self):
        """Correct spectra for stray light - spectra are assumed to be dark-corrected prior to running this function"""
        # Set the range of stray pixels
        stray_range = np.arange(self._start_stray_pix, self._end_stray_pix + 1)

        # Correct clear and plume spectra (assumed to already be dark subtracted)
        if not self.stray_corrected_clear:
            self.clear_spec_corr = self.clear_spec_corr - np.mean(self.clear_spec_corr[stray_range])
            self.clear_spec_corr[self.clear_spec_corr < 0] = 0
            self.stray_corrected_clear = True

        if not self.stray_corrected_plume:
            self.plume_spec_corr = self.plume_spec_corr - np.mean(self.plume_spec_corr[stray_range])
            self.plume_spec_corr[self.plume_spec_corr < 0] = 0
            self.stray_corrected_plume = True

    def load_ref_spec(self, pathname, species):
        """Load raw reference spectrum"""
        self.ref_spec[species] = np.loadtxt(pathname)

        # Remove un-needed data above 400 nm, which may slow down processing
        idxs = np.where(self.ref_spec[species][:, 0] > 400)
        self.ref_spec[species] = self.ref_spec[species][:np.min(idxs), :]

        # Assume we have loaded a new spectrum, so set this to False - ILS has not been convolved yet
        self.ref_convolved = False


    def get_ref_spectrum(self):
        """Load in reference spectrum"""
        self.wavelengths = None  # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        #
        # --------------------------------

    def conv_ref_spec(self):
        """
        Convolves reference spectrum with instument line shape (ILS)
        after first interpolating to spectrometer wavelengths
        """
        if self.wavelengths is None:
            print('No wavelength data to perform convolution')
            return

        # Need an odd sized array for convolution, so if even we omit the last pixel
        if self.ILS.size % 2 == 0:
            self.ILS = self.ILS[:-1]

        # Loop through all reference species we have loaded and resample their data to the spectrometers wavelengths
        for f in self.ref_spec_used:
            self.ref_spec_interp[f] = np.interp(self.wavelengths, self.ref_spec[f][:, 0], self.ref_spec[f][:, 1])
            self.ref_spec_conv[f] = convolve(self.ref_spec_interp[f], self.ILS)
            self.ref_spec_ppmm[f] = self.ref_spec_conv[f] * self.ppmm_conversion

        # Update bool as we have now performed this process
        self.ref_convolved = True

    def load_calibration_spectrum(self, pathname):
        """Load Calibation image for spectrometer"""
        pass

    def set_fit_windows(self):
        """Define fitting window for DOAS procedure
        If wavelength domain is used, first convert this to pixel space"""
        if self.wave_fit:
            if self.wavelengths is None:
                print('Error, first run get_ref_spectrum() to define wavelengths vector')
                return

        self.fit_window = np.arange(self._start_fit_pix, self._end_fit_pix)  # Fitting window (in Pixel space)
        self.fit_window_ref = self.fit_window - self.shift

    def stretch_spectrum(self, ref_key):
        """Stretch/squeeze reference spectrum to improve fit"""
        if self.stretch == 0:
            # If no stretch, extract reference spectrum using fit window and return
            return self.ref_spec_filter[ref_key][self.fit_window_ref]
        else:
            stretch_inc = (self.stretch * self.stretch_adjuster) / self.stretch_resample  # Stretch increment for resampled spectrum

            if self.stretch < 0:
                # Generate the fit window for extraction
                # Must be larger than fit window as a squeeze requires pulling in data outside of the fit window
                if self.fit_window_ref[-1] < (len(self.wavelengths) - 50):
                    extract_window = np.arange(self.fit_window_ref[0], self.fit_window_ref[-1] + 50)
                else:
                    extract_window = np.arange(self.fit_window_ref[0], len(self.wavelengths))
            else:
                extract_window = self.fit_window_ref

            wavelengths = self.wavelengths[extract_window]
            values = self.ref_spec_filter[ref_key][extract_window]

            # Generate new arrays with 'stretch_resample' more data points
            wavelengths_resampled = np.linspace(wavelengths[0], wavelengths[-1], len(wavelengths)*self.stretch_resample)
            values_resample = np.interp(wavelengths_resampled, wavelengths, values)

            num_pts = len(wavelengths_resampled)
            wavelengths_stretch = np.zeros(num_pts)

            # Stretch wavelengths
            for i in range(num_pts):
                wavelengths_stretch[i] = wavelengths_resampled[i] + (i * stretch_inc)

            values_stretch = np.interp(self.wavelengths[self.fit_window_ref], wavelengths_stretch, values_resample)

            return values_stretch

    def load_spec(self):
        """Load spectrum"""
        pass

    def load_dark(self):
        """Load drk images -> co-add to generate single dark image"""
        pass

    def save_dark(self, filename):
        """Save dark spectrum"""
        if self.wavelengths is None or self.dark_spec is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.dark_spec):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.dark_spec]),
                   header='Dark spectrum\nWavelength [nm]\tIntensity [DN]')

    def save_clear_raw(self, filename):
        """Save clear spectrum"""
        if self.wavelengths is None or self.clear_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.clear_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.clear_spec_raw]),
                   header='Raw clear (Fraunhofer) spectrum\n'
                          '-Not dark-corrected\nWavelength [nm]\tIntensity [DN]')

    def save_plume_raw(self, filename):
        if self.wavelengths is None or self.plume_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.plume_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.plume_spec_raw]),
                   header='Raw in-plume spectrum\n'
                          '-Not dark-corrected\nWavelength [nm]\tIntensity [DN]')

    @staticmethod
    def doas_fit(ref_spec, *fit_params):
        """
        DOAS fit function

        Parameters
        -----------
        ref_spec : array-like object
            Contains all (k) reference spectra to be fitted in a (k,N) array.
            N is the size of the spectrum grid
        """
        # Unpack fit scalars into column vector
        scalars = np.transpose(np.array([list(fit_params)]))

        # Sum reference spectra to create total absorbance for comparison with measured absorbance spectrum
        total_absorbance = np.sum(ref_spec * scalars, axis=0)

        return total_absorbance

    def poly_doas(self):
        """
        Performs main processing in polynomial fitting DOAS retrieval
        NOT COMPLETE, DO NOT USE!!!
        """
        self.wavelengths_cut = self.wavelengths[self.fit_window]  # Extract wavelengths (used in plotting)

        with np.errstate(divide='ignore'):
            self.abs_spec = np.log(np.divide(self.clear_spec_corr, self.plume_spec_corr))  # Calculate absorbance

        self.ref_spec_cut['SO2'] = self.ref_spec_ppmm[self.ref_spec_types['SO2']][self.fit_window_ref]
        self.abs_spec_cut = self.abs_spec[self.fit_window]

        idx = 0
        for i in self.vals_ca:
            ref_spec_fit = self.ref_spec_scaled['SO2'] * i  # Our iterative guess at the SO2 column density
            residual = self.abs_spec_cut - ref_spec_fit  # Calculate resultant residual from spectrum fitting
            poly_fit = np.polyfit(self.fit_window, residual, self.poly_order)  # Fit polynomial to residual
            poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window

            self.mse_vals[idx] = np.mean(np.power(residual - poly_vals, 2))  # Calculate MSE of fit

            idx += 1

        self.min_idx = np.argmin(self.mse_vals)
        self.column_density['SO2'] = self.vals_ca[self.min_idx]

    def fltr_doas(self):
        """
        Performs main retrieval in digital filtering DOAS retrieval
        """
        self.wavelengths_cut = self.wavelengths[self.fit_window]    # Extract wavelengths (used in plotting)

        # Calculate absorbance
        with np.errstate(divide='ignore', invalid='ignore'):
            self.abs_spec = np.log(np.divide(self.clear_spec_corr, self.plume_spec_corr))

        # Remove nans and infs (filter function doesn't work if they are included)
        self.abs_spec[np.isinf(self.abs_spec)] = np.nan
        self.abs_spec[np.isnan(self.abs_spec)] = 0

        self.abs_spec_filt = signal.lfilter(self.filt_B, self.filt_A, self.abs_spec)    # Filter absorbance spectrum
        self.abs_spec_cut = self.abs_spec_filt[self.fit_window]                         # Extract fit window

        # Loop through reference spectra and prepare them for processing
        for spec in self.ref_spec_used:

            # Filter reference spectrum and extract fit window
            self.ref_spec_filter[spec] = signal.lfilter(self.filt_B, self.filt_A, self.ref_spec_ppmm[spec])

            # Stretch spectrum
            self.ref_spec_cut[spec] = self.stretch_spectrum(spec)

        # self.ref_spec_cut['SO2'] = self.ref_spec_filter['SO2'][self.fit_window_ref]

        # ------------------------------------------------------------------------------------------
        # Scipy.optimize.curve_fit least squares fitting
        # ------------------------------------------------------------------------------------------
        # Pack all requested reference spectra into an array for curve fitting
        ref_spectra_packed = np.empty((len(self.ref_spec_used), len(self.abs_spec_cut)))
        i = 0
        for spec in self.ref_spec_used:
            ref_spectra_packed[i, :] = self.ref_spec_cut[spec]
            i += 1

        # Run fit
        column_densities, pcov = curve_fit(self.doas_fit, ref_spectra_packed, self.abs_spec_cut,
                                           p0=np.ones(ref_spectra_packed.shape[0]))
        self.std_err = round(np.sqrt(np.diag(pcov))[0], 1)

        # Loop through species to unpack results
        i = 0
        self.ref_spec_fit['Total'] = np.zeros(len(self.abs_spec_cut))
        for spec in self.ref_spec_used:
            # Unpack column densities to dictionary
            self.column_density[spec] = int(round(column_densities[i]))

            # Generate scaled reference spectrum
            self.ref_spec_fit[spec] = self.ref_spec_cut[spec] * self.column_density[spec]

            # Add scaled reference spectrum to total ref spec
            self.ref_spec_fit['Total'] += self.ref_spec_fit[spec]

            i += 1

    def gen_abs_specs(self):
        """
        Generate absorbance spectra for individual species and match to measured absorbance
        :return:
        """
        # Iterate through reference spectra
        for spec in self.ref_spec_used:

            # Retrieve all ref spectra other than that held in spec
            other_species = [f for f in self.ref_spec_used if f is not spec]

            # Subtract contributions of these other spectra to absorption
            self.abs_spec_species[spec] = self.abs_spec_cut
            for i in other_species:
                self.abs_spec_species[spec] = self.abs_spec_species[spec] - self.ref_spec_fit[i]

        # Setup total absorbance spectrum as 'Total' key
        self.abs_spec_species['Total'] = self.abs_spec_cut

        # Calculate final residual by subtracting all absorbing species
        self.abs_spec_species['residual'] = self.abs_spec_cut
        for spec in self.ref_spec_used:
            self.abs_spec_species['residual'] = self.abs_spec_species['residual'] - self.ref_spec_fit[spec]


    def poly_plot_gen(self):
        """Generate arrays to be plotted -> residual, fitted spectrum"""
        self.ref_spec_fit['SO2'] = self.ref_spec_cut['SO2'] * self.column_density['SO2']
        self.residual = self.abs_spec_cut - self.ref_spec_fit['SO2']
        poly_fit = np.polyfit(self.fit_window, self.residual, self.poly_order)  # Fit polynomial to residual
        self.poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window
        self.best_fit = self.ref_spec_fit['SO2'] + self.poly_vals  # Generate best fit absorbance spectrum

        # MAKE PLOT
        plt.figure()
        abs_plt, = plt.plot(self.abs_spec_cut, label='Absorbance spectrum')
        ref_plt, = plt.plot(self.ref_spec_fit, label='Reference spectrum * CA')
        res_plt, = plt.plot(self.residual, label='Residual')
        poly_plt, = plt.plot(self.poly_vals, label='Polynomial fit')
        best_plt, = plt.plot(self.best_fit, label='Best fit')
        plt.xlabel('Pixel')
        plt.ylabel('Absorbance')
        plt.legend(handles=[abs_plt, ref_plt, res_plt, poly_plt, best_plt])
        plt.show()

    def process_doas(self):
        """Handles the order of DOAS processing"""
        # Check we have all of the correct spectra to perform processing
        if self.clear_spec_raw is None or self.plume_spec_raw is None or self.wavelengths is None:
            raise SpectraError('Require clear and plume spectra for DOAS processing')

        if self.ref_spec_types[0] not in self.ref_spec.keys():
            raise SpectraError('No SO2 reference spectrum present for processing')

        # If the stray region hasn't been converted to pixel space we just set it to itself to implement pixel mapping
        if self._start_stray_pix is None or self._end_stray_pix is None:
            self.start_stray_wave = self.start_stray_wave
            self.end_stray_wave = self.end_stray_wave

        # Same for fit window
        if self._start_fit_pix is None or self._end_fit_pix is None:
            self.start_fit_wave = self.start_fit_wave
            self.end_fit_wave = self.end_fit_wave

        if not self.dark_corrected_clear or not self.dark_corrected_plume:
            if self.dark_spec is None:
                print('Warning! No dark spectrum present, processing without dark subtraction')

                # Set raw spectra to the corrected spectra, ignoring that they have not been dark corrected
                self.clear_spec_corr = self.clear_spec_raw
                self.plume_spec_corr = self.plume_spec_raw
            else:
                self.dark_corr_spectra()

        # Correct spectra for stray light
        if not self.stray_corrected_clear or not self.stray_corrected_plume:
            self.stray_corr_spectra()

        # Set fitting windows for acquired and reference spectra
        self.set_fit_windows()

        # Convolve reference spectrum with the instrument lineshape
        if not self.ref_convolved:
            self.conv_ref_spec()

        # Run processing
        self.fltr_doas()

        # Generate absorbance spectra for individual species
        self.gen_abs_specs()

        # Set flag defining that data has been fully processed
        self.processed_data = True


class SpectraError(Exception):
    """
    Error raised if correct spectra aren't present for processing
    """
    pass


class SpectrometerCal:
    """Class to calibrate spectrometer"""
    def __init__(self):
        pass


class ScanProcess:
    """
    Class to control processing of DOAS scan data
    """
    def __init__(self):
        self.plume_distance = None  # Distance to plume [m]
        self.plume_speed = None     # Plume speed [m/s]
        self.scan_sep = None       # Distance between two scan points in the plume [m]
        self.scan_incr = ScanProperties().scan_incr     # Scanner angle increment

        self.ppm2kg = 2.663 * 1.0e-6    # Conversion factor to get ppm.m in kg/m2

        self.scan_angles = np.array([])
        self.column_densities = np.array([])

        self.SO2_flux = None

    def clear_data(self):
        """Initialises new arrays"""
        self.scan_angles = np.array([])
        self.column_densities = np.array([])

    def add_data(self, scan_angle, column_density):
        """Adds data"""
        self.scan_angles = np.append(self.scan_angles, scan_angle)
        self.column_densities = np.append(self.column_densities, column_density)

    def __calc_scan_sep__(self):
        """Calculates separation between 2 scan points in the plume"""
        scan_step = np.deg2rad(self.scan_angles[1] - self.scan_angles[0])
        self.scan_sep = 2 * np.tan(scan_step/2) * self.plume_distance

    def calc_emission_rate(self):
        """Calculates emission rate from current data"""
        self.__calc_scan_sep__()

        # Convert column densities to kg/m2
        cd_conv = self.column_densities * self.ppm2kg

        # Find total mass of SO2 across scan profile kg/m
        SO2_mass = integrate.trapz(cd_conv, dx=self.scan_sep)

        # Calculate emssion rate (flux) kg/s
        self.SO2_flux = SO2_mass * self.plume_speed

    @property
    def flux_tons(self):
        """Convert SO2 flux to t/day"""
        if self.SO2_flux:
            return (self.SO2_flux/1000) * 60 * 60 * 24
        else:
            return None


class EmissionSeries:
    """
    Holds emission rate timeseries data
    """
    def __init__(self):
        self.datetime_fmt ='%Y-%m-%d\t%H%M%S'
        self.filename = 'emission_series.txt'
        self.times = []     # List of times (not numpy array as times are held as datetime objects)
        self.emission_rates = np.array([])

    def clear_data(self):
        """Initialises new arrays"""
        self.times = []
        self.emission_rates = np.array([])

    def add_data(self, dat_time, emission_rate):
        """Adds data"""
        self.times.append(dat_time)
        self.emission_rates = np.append(self.emission_rates, emission_rate)

    @property
    def matplotlib_times(self):
        return date2num(self.times)

    @property
    def str_times(self):
        return [d.strftime(self.datetime_fmt) for d in self.times]

    @property
    def emission_tons(self):
        """Returns emission rate in t/day"""
        if len(self.emission_rates) > 0:
            return (self.emission_rates/1000) * 60 * 60 * 24
        else:
            return None

    @property
    def mean_emission(self):
        if len(self.emission_rates) > 0:
            return np.mean(self.emission_rates)
        else:
            return None

    @property
    def max_emission(self):
        if len(self.emission_rates) > 0:
            return np.max(self.emission_rates)
        else:
            return None

    def save_series(self, pathname, plume_distance, plume_speed):
        """
        Saves times series data
        :param pathname:    str     Path to save file directory
        :return:
        """
        # Create full path
        pathname = os.path.join(pathname, self.filename)

        # Sometimes encounter an exception when trying to save, issue with permission which seems random, so keep trying
        # to save until we are successful
        while True:
            try:
                np.savetxt(pathname, np.transpose([self.str_times, self.emission_rates]),
                           header='plume_distance={}\n'
                                  'plume_speed={}\n'
                                  'Date\tTime\tEmission rate [kg/s]'.format(plume_distance, plume_speed))
                break
            except Exception as e:
                pass

if __name__ == "__main__":
    doas_process = DOASWorker(2)
    doas_process.get_ref_spectrum()
    doas_process.set_fit_windows()