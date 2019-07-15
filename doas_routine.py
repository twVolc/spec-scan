# Main Subroutine which processes images according to the DOAS retrieval method.

# THIS IS TOO MESSY, I NEED TO SEPARATE THE GUI FROM THE DOAS PROCESSING!

import numpy as np
from scipy import signal
import os
import glob
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from astropy.convolution import convolve
import scipy.integrate as integrate

class DOASWorker:
    """Class to control DOAS processing
    General order of play for processing:
    Initiate class,
    get_ref_spectrum()
    set_fit_window()
    shift_spectrum()"""
    def __init__(self, routine=2):
        self.routine = routine  # Defines routine to be used, either (1) Polynomial or (2) Digital Filtering

        # ======================================================================================================================
        # Initial Definitions
        # ======================================================================================================================
        self.ppmm_conversion = 2.7e15   # convert absorption cross-section in cm2/molecule to ppm.m (MAY NEED TO CHANGE THIS TO A DICTIONARY AS THE CONVERSION MAY DIFFER FOR EACH SPECIES?)

        self.shift = 0                  # Shift of spectrum in number of pixels
        self.stretch = 0                # Stretch of spectrum
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
        self.abs_spec = None
        self.abs_spec_cut = None
        self.abs_spec_filt = None
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

        # Assume we have loaded a new spectrum, so set this to False - ILS has not been convolved yet
        self.ref_convolved = False


    def get_ref_spectrum(self):
        """Load in reference spectrum"""
        self.wavelengths = None  # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        #
        # --------------------------------

    def conv_ref_spec(self):
        """Convolves reference spectrum with instument line shape (ILS)
        after first interpolating to spectrometer wavelengths"""
        if self.wavelengths is None:
            print('No wavelength data to perform convolution')
            return

        species = [f for f in self.ref_spec_types if f in self.ref_spec.keys()]

        # Need an odd sized array for convolution, so if even we omit the last pixel
        if self.ILS.size % 2 == 0:
            self.ILS = self.ILS[:-1]

        # Loop through all reference species we have loaded and resample their data to the spectrometers wavelengths
        for f in species:
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
        self.column_amount = self.vals_ca[self.min_idx]

    def fltr_doas(self):
        """Performs main retrieval in digital filtering DOAS retrieval"""
        self.wavelengths_cut = self.wavelengths[self.fit_window]    # Extract wavelengths (used in plotting)

        with np.errstate(divide='ignore', invalid='ignore'):
            self.abs_spec = np.log(np.divide(self.clear_spec_corr, self.plume_spec_corr))   # Calculate absorbance

        # Remove nans and infs (filter function doesn't work if they are included)
        self.abs_spec[np.isinf(self.abs_spec)] = np.nan
        self.abs_spec[np.isnan(self.abs_spec)] = 0

        self.abs_spec_filt = signal.lfilter(self.filt_B, self.filt_A, self.abs_spec)    # Filter absorbance spectrum
        self.abs_spec_cut = self.abs_spec_filt[self.fit_window]                         # Extract fit window

        # Filter reference spectrum and extract fit window
        self.ref_spec_filter['SO2'] = signal.lfilter(self.filt_B, self.filt_A, self.ref_spec_ppmm['SO2'])
        # PUT STRETCH HERE
        self.ref_spec_cut[self.ref_spec_types[0]] = self.stretch_spectrum(self.ref_spec_types[0])
        # self.ref_spec_cut['SO2'] = self.ref_spec_filter['SO2'][self.fit_window_ref]

        # ------------------------------------------------------------------------------------------
        # Attempting faster iterative process by not using every vals_ca value initially
        idx = 0
        for i in self.vals_ca_cut:
            ref_spec_fit = self.ref_spec_cut['SO2'] * i

            self.mse_vals_cut[idx] = np.mean(np.power(self.abs_spec_cut - ref_spec_fit, 2))  # Calculate MSE of fit
            idx += 1

        # Find best fit, then hone in on that region to iterate through at a step of 1 ppm.m
        min_idx_1 = np.argmin(self.mse_vals_cut)

        if min_idx_1 == 0:
            vals_ca_new = self.vals_ca[:self.vals_ca_cut_idxs[min_idx_1 + 2]]
        elif min_idx_1 == len(self.mse_vals_cut) - 1:
            vals_ca_new = self.vals_ca[self.vals_ca_cut_idxs[min_idx_1 - 2]:]
        else:
            vals_ca_new = self.vals_ca[self.vals_ca_cut_idxs[min_idx_1 - 2]:self.vals_ca_cut_idxs[min_idx_1 + 2]]

        mse_vals_new = np.zeros(len(vals_ca_new))
        idx = 0
        for i in vals_ca_new:
            ref_spec_fit = self.ref_spec_cut['SO2'] * i

            mse_vals_new[idx] = np.mean(np.power(self.abs_spec_cut - ref_spec_fit, 2))  # Calculate MSE of fit
            idx += 1
        min_idx_2 = np.argmin(mse_vals_new)
        self.column_amount = vals_ca_new[min_idx_2]


        # Standard slow iterative process
        # idx = 0
        # for i in self.vals_ca:
        #     ref_spec_fit = self.ref_spec_cut['SO2'] * i
        #
        #     self.mse_vals[idx] = np.mean(np.power(self.abs_spec_cut - ref_spec_fit, 2))  # Calculate MSE of fit
        #     idx += 1

        # Determine column amount from best fit
        # self.min_idx = np.argmin(self.mse_vals)
        # self.column_amount = self.vals_ca[self.min_idx]

        # Generate scaled reference spectrum
        self.ref_spec_fit['SO2'] = self.ref_spec_cut['SO2'] * self.column_amount


    def poly_plot_gen(self):
        """Generate arrays to be plotted -> residual, fitted spectrum"""
        self.ref_spec_fit['SO2'] = self.ref_spec_cut['SO2'] * self.column_amount
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


if __name__ == "__main__":
    doas_process = DOASWorker(2)
    doas_process.get_ref_spectrum()
    doas_process.set_fit_windows()