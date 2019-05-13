# Main Subroutine which processes images according to the DOAS retrieval method.

# THIS IS TOO MESSY, I NEED TO SEPARATE THE GUI FROM THE DOAS PROCESSING!

import numpy as np
from scipy import signal
import os
import glob
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

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
        self.shift = 0              # Shift of spectrum in number of pixels
        self.stretch = 0            # Stretch of spectrum
        self._start_stray_pix = None  # Pixel space stray light window definitions
        self._end_stray_pix = None
        self._start_stray_wave = 280  # Wavelength space stray light window definitions
        self._end_stray_wave = 290
        self._start_fit_pix = None  # Pixel space fitting window definitions
        self._end_fit_pix = None
        self._start_fit_wave = 305  # Wavelength space fitting window definitions
        self._end_fit_wave = 320
        self.fit_window = None      # Fitting window, determined by set_fit_window()
        self.fit_window_ref = None  # Placeholder for shifted fitting window for the reference spectrum
        self.wave_fit = True        # If True, wavelength parameters are used to define fitting window

        self.wavelengths = None     # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        self.dark_spec = None       # Dark spectrum
        self.clear_spec_raw = None  # Clear (fraunhofer) spectrum - not dark corrected
        self.plume_spec_raw = None  # In-plume spectrum (main one which is used for calculation of SO2
        self.ref_spec = dict()      # Create empty dictionary for holding reference spectra
        self.ref_spec_types = ['SO2', 'O3', 'ring'] # List of reference spectra types accepted/expected

        self.poly_order = 2  # Order of polynomial used to fit residual
        (self.filt_B, self.filt_A) = signal.butter(10, 0.065, btype='highpass')

        self.start_ca = -1000  # Starting column amount for iterations
        self.end_ca = 10000  # Ending column amount for iterations
        self.vals_ca = np.arange(self.start_ca, self.end_ca+1)  # Array of column amounts to be iterated over
        self.mse_vals = np.zeros(len(self.vals_ca))  # Array to hold mse values

        self.filetypes = dict(defaultextension='.png', filetypes=[('PNG', '*.png')])

        # ----------------------------------------------------------------------------------
        # We need to make sure that all images are dark subtracted before final processing
        # Also make sure that we don't dark subtract more than once!
        self.have_dark = False  # Used to define if a dark image is loaded.
        self.cal_dark_corr = False  # Tells us if the calibration image has been dark subtracted
        self.clear_dark_corr = False  # Tells us if the clear image has been dark subtracted
        self.plume_dark_corr = False  # Tells us if the plume image has been dark subtracted

        # ==============================================================================================================

        # # --------------------------------------------------------------------------------------------------------------
        # # GENERATE CLEAR SPECTRUM AND LOAD DARK IMAGE
        # (self.img_clear, self.img_size_x, self.img_size_y) = self.load_spec()  # Clear image (I0)
        # self.img_dark = self.load_dark()  # Dark Image
        # self.img_clear = self.img_clear - self.img_dark  # Dark subtract clear image
        # # --------------------------------------------------------------------------------------------------------------

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



    def load_reference_spectrum(self, pathname, species):
        """Load raw reference spectrum"""
        self.ref_spec[species] = np.loadtxt(pathname)


    def get_ref_spectrum(self):
        """Load in reference spectrum"""
        self.wavelengths = None  # Placeholder for wavelengths attribute which contains all wavelengths of spectra
        #
        # --------------------------------

    def load_calibration_spectrum(self, pathname):
        """Load Calibation image for spectrometer"""
        pass

    def set_fit_window(self):
        """Define fitting window for DOAS procedure
        If wavelength domain is used, first convert this to pixel space"""
        if self.wave_fit:
            if self.wavelengths is None:
                print('Error, first run get_ref_spectrum() to define wavelengths vector')
                return

            ## THis shouldn't be necessary now as I set these during setting the wavelength fit start/end
            # wave_dif_start = self.wavelengths - self.start_fit_wave
            # wave_dif_end = self.wavelengths - self.end_fit_wave
            # self._start_fit_pix = wave_dif_start.index(np.amin(wave_dif_start))  # Find the index which represents the wavelengths closest to the defined starting wavelength for the fit
            # self._end_fit_pix = wave_dif_end.index(np.amin(wave_dif_end))  # As above, but for ending wavelength

        self.fit_window = np.arange(self._start_fit_pix, self._end_fit_pix)  # Fitting window (in Pixel space)

    def shift_spectrum(self):
        """Shift fitting window for reference spectrum"""
        self.fit_window_ref = self.fit_window - self.shift  # Shifting the fitting window indices for the ref spectrum

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
        """Save dark spectrum"""
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


    def poly_DOAS(self):
        """Performs main processing in polynomial fitting DOAS retrieval"""

        self.abs_spec = np.log(np.divide(self.clear_spec_raw, self.plume_spec))  # Calculate absorbance
        self.ref_spec_cut = self.ref_spec[self.fit_window_ref]
        self.abs_spec_cut = self.abs_spec[self.fit_window]

        idx = 0
        for i in self.vals_ca:
            ref_spec_fit = self.ref_spec_cut * i  # Our iterative guess at the SO2 column density
            residual = self.abs_spec_cut - ref_spec_fit  # Calculate resultant residual from spectrum fitting
            poly_fit = np.polyfit(self.fit_window, residual, self.poly_order)  # Fit polynomial to residual
            poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window

            self.mse_vals[idx] = np.mean(np.power(residual - poly_vals, 2))  # Calculate MSE of fit

            idx += 1

        self.min_idx = np.argmin(self.mse_vals)
        self.column_amount = self.vals_ca[self.min_idx]

    def fltr_DOAS(self):
        """Performs main retrieval in digital filtering DOAS retrieval"""
        self.abs_spec = np.log(np.divide(self.clear_spec_raw, self.plume_spec))  # Calculate absorbance
        self.abs_spec_filt = signal.lfilter(self.filt_B, self.filt_A, self.abs_spec)  # Filter absorbance spectrum

        self.ref_spec_cut = self.ref_spec[self.fit_window_ref]
        self.abs_spec_cut = self.abs_spec[self.fit_window]

        idx = 1
        for i in self.vals_ca:
            ref_spec_fit = self.ref_spec_cut * i

            self.mse_vals[idx] = np.mean(np.power(self.abs_spec_cut, ref_spec_fit, 2))  # Calculate MSE of fit
            idx += 1

        self.min_idx = np.argmin(self.mse_vals)
        self.column_amount = self.vals_ca[self.min_idx]



    def poly_plot_gen(self):
        """Generate arrays to be plotted -> residual, fitted spectrum"""
        self.ref_spec_fit = self.ref_spec_cut * self.column_amount
        self.residual = self.abs_spec_cut - self.ref_spec_fit
        poly_fit = np.polyfit(self.fit_window, self.residual, self.poly_order)  # Fit polynomial to residual
        self.poly_vals = np.polyval(poly_fit, self.fit_window)  # Generate polynomial values for fitting window
        self.best_fit = self.ref_spec_fit + self.poly_vals  # Generate best fit absorbance spectrum

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


    def process_DOAS(self):
        """Handles the order of DOAS processing"""
        # Check we have all of the correct spectra to perform processing
        if self.clear_spec_raw is None or self.plume_spec_raw is None or self.wavelengths:
            raise SpectraError('Require clear and plume spectra for DOAS processing')

        if self.ref_spec_types[0] not in self.ref_spec.keys():
            raise SpectraError('No SO2 reference spectrum present for processing')

        if self.dark_spec is None:
            print('Warning! No dark spectrum present, processing without dark subtraction')







class SpectraError(Exception):
    """
    Error raised if correct spectra aren't present for processing
    """
    pass


class SpectrometerCal:
    """Class to calibrate spectrometer"""
    def __init__(self):
        pass


if __name__ == "__main__":
    doas_process = DOASWorker(2)
    doas_process.get_ref_spectrum()
    doas_process.set_fit_window()
    doas_process.shift_spectrum()