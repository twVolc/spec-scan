# -*- coding: utf-8 -*-

"""
Main controller classes for the PiCam and OO Flame spectrometer
"""

import warnings
import queue
import multiprocessing.managers
import io
import os
import time
import datetime
import numpy as np
import threading
import serial

from .setupclasses import SpecSpecs
from .utils import format_time

try:
    import seabreeze.spectrometers as sb
except ModuleNotFoundError:
    warnings.warn('Working on machine without seabreeze, functionality of some classes will be lost')


class Spectrometer(SpecSpecs):
    """Main class for spectrometer control

    subclass of :class: SpecSpecs
    """
    def __init__(self, filename=None):
        super().__init__(filename)

        self.capture_q = queue.Queue()      # Queue for requesting spectra
        self.spec_q = queue.Queue()         # Queue to put spectra in for access elsewhere
        self.capture_thread = None          # Thread for interactive capture

        # Discover spectrometer devices
        self.devices = None  # List of detected spectrometers
        self.spec = None  # Holds spectrometer for interfacing via seabreeze
        self.find_device()

        # Set integration time (ALL IN MICROSECONDS)
        self._int_limit_lower = 1000  # Lower integration time limit
        self._int_limit_upper = 20000000  # Upper integration time limit
        self._int_time = None  # Integration time attribute
        self.int_time = self.start_int_time

        self.min_coadd = 1
        self.max_coadd = 100
        self._coadd = None  # Controls coadding of spectra
        self.coadd = self.start_coadd

        self.continuous_capture = False     # Bool set to true when camera is in continuous capture mode

    def find_device(self):
        """Function to search for devices"""
        try:
            self.devices = sb.list_devices()
            self.spec = sb.Spectrometer(self.devices[0])
            self.spec.trigger_mode(0)

            # If we have a spectrometer we then retrieve its wavelength calibration and store it as an attribute
            self.get_wavelengths()

        except IndexError:
            self.devices = None
            self.spec = None
            raise SpectrometerConnectionError('No spectrometer found')

    def _q_check(self, q, q_type='capt'):
        """Checks type of queue object and returns queue (ret_q). Sets queue to default queue if none is provided"""
        if isinstance(q, multiprocessing.managers.BaseProxy):
            print('Using multiprocessing queue')
            ret_q = q
        elif isinstance(q, queue.Queue):
            print('Using Queue queue')
            ret_q = q
        else:
            print('Unrecognized queue object, reverting to default')
            if q_type == 'capt':
                ret_q = self.capture_q
            elif q_type == 'spec':
                ret_q = self.spec_q
            else:
                ret_q = queue.Queue()

        return ret_q

    @property
    def int_time(self):
        return self._int_time / 1000  # Return time in milliseconds

    @int_time.setter
    def int_time(self, int_time):
        """Set integration time

        Parameters
        ----------
        int_time: int
            Integration time for spectrometer, provided in milliseconds
        """
        # Adjust to work in microseconds (class takes time in milliseconds) and ensure we have an <int>
        int_time = int(int_time * 1000)

        # Check requested integration time is acceptable
        if int_time < self._int_limit_lower:
            raise ValueError('Integration time below %i us is not possible' % self._int_limit_lower)
        elif int_time > self._int_limit_upper:
            raise ValueError('Integration time above %i us is not possible' % self._int_limit_upper)

        self._int_time = int_time

        # Adjust _int_time_idx to reflect the closest integration time to the current int_time
        self._int_time_idx = np.argmin(np.abs(self.int_list - self.int_time))

        # Set spectrometer integration time
        self.spec.integration_time_micros(int_time)

    @property
    def int_time_idx(self):
        return self._int_time_idx

    @int_time_idx.setter
    def int_time_idx(self, value):
        """Update integration time to value in int_list defined by int_time_idx when int_time_idx is changed"""
        self._int_time_idx = value
        self.int_time = self.int_list[self.int_time_idx]

    @property
    def coadd(self):
        return self._coadd

    @coadd.setter
    def coadd(self, coadd):
        """Set coadding property"""
        if coadd < self.min_coadd:
            coadd = self.min_coadd
        elif coadd > self.max_coadd:
            coadd = self.max_coadd
        self._coadd = int(coadd)

    def generate_filename(self, time_str, spec_type):
        """Generates the spectrum filename

        Parameters
        ----------
        time_str: str
            Time string containing date and time
        """
        return time_str + '_' + self.file_ss.format(int(self.int_time)) + '_' \
               + str(self.coadd) + 'coadd_' + spec_type + self.file_ext

    def get_spec(self):
        """Acquire spectrum from spectrometer"""
        # Set array for coadding spectra
        coadded_spectrum = np.zeros(len(self.wavelengths))

        # Loop through number of coadds
        for i in range(self.coadd):
            coadded_spectrum += self.spec.intensities()

        # Correct for number of coadds to result in a spectrum with correct digital numbers for bit-depth of device
        coadded_spectrum /= self.coadd
        self.spectrum = coadded_spectrum

    def get_spec_now(self):
        """Immediately acquire spectrum from spectrometer - does not discard first spectrum (probably never used)"""
        self.spectrum = self.spec.intensities()

    def get_wavelengths(self):
        """Returns wavelengths"""
        self.wavelengths = self.spec.wavelengths()

    def extract_subspec(self, wavelengths):
        """Extract and return wavelengths and spectrum data for subsection of spectrum defined by wavelengths

        Parameters
        ----------
        wavelengths: list, tuple

        Returns
        -------
        wavelengths: list
            wavelengths of spectrometer extracted between range requested
        spectrum: list
            intensities from spectrum extracted between requested range
        """
        # Check wavelengths have been provided correctly
        if len(wavelengths) != 2:
            raise ValueError('Expected list or tuple of length 2')

        # Determine indices of arrays where wavelengths are closest to requested extraction wavelengths
        min_idx = np.argmin(np.abs(wavelengths[0] - self.wavelengths))
        max_idx = np.argmin(np.abs(wavelengths[1] - self.wavelengths))

        return self.wavelengths[min_idx:max_idx+1], self.spectrum[min_idx:max_idx+1]

    def check_saturation(self):
        """Check spectrum saturation
        return -1: if saturation exceeds the maximum allowed
        return 1:  if saturation is below minimum allowed
        return 0:  otherwise
        """
        # Extract spectrum in specific wavelength range to be checked
        wavelengths, spectrum = self.extract_subspec(self.saturation_range)

        if np.amax(spectrum) / self._max_DN > self.max_saturation:
            return -1
        elif np.amax(spectrum) / self._max_DN < self.min_saturation:
            return 1
        else:
            return 0

    def interactive_capture(self, spec_q=None, capt_q=None):
        """Public access thread starter for _interactive_capture()"""
        self.capture_thread = threading.Thread(target=self._interactive_capture, args=(spec_q, capt_q,))
        self.capture_thread.daemon = True
        self.capture_thread.start()

    def _interactive_capture(self, spec_q=None, capt_q=None):
        """Interactive capturing by requesting captures through capt_q

        Parameters
        ---------
        spec_q: Queue-like object
            Spectra are passed to this queue once captured
        capt_q: Queue-like object
            Capture commands are passed to this object using its put() method
        """

        # Setup queue
        capt_q = self._q_check(capt_q, q_type='capt')
        spec_q = self._q_check(spec_q, q_type='spec')

        while True:

            # Wait for imaging command (expecting a dictionary containing information for acquisition)
            command = capt_q.get(block=True)

            if 'exit' in command:
                # return if commanded to exit
                if command['exit']:
                    return

            if 'int_time' in command:
                # Set shutter speed
                self.int_time = command['int_time']

            # Start a continous capture if requested
            if 'start_cont' in command:
                if command['start_cont']:
                    # If we have been provided with a queue for images we pass this to capture_sequence()
                    if 'spec_q' in command:
                        self.capture_sequence(spec_q=command['spec_q'], capt_q=capt_q)
                    else:
                        self.capture_sequence(spec_q=spec_q, capt_q=capt_q)

            # Instigate capture of all dark images
            elif 'dark_seq' in command:
                if command['dark_seq']:
                    self.capture_darks()

            # If continuous capture is not requested we check if any single image is requested
            else:
                if 'type' in command:
                    # If a sequence isn't requested we take one typical image
                    if command['type'] in self.file_spec_type:

                        # Get time and format
                        time_str = format_time(datetime.datetime.now(), self.file_datestr)

                        # Capture spectrum
                        self.get_spec()

                        # Generate filename
                        filename = self.generate_filename(time_str, command['type'])

                        # Put filename and spectrum in queue
                        spec_q.put([filename, self.spectrum])

    def capture_sequence(self, spec_q=None, capt_q=None):
        """Captures sequence of spectra

        Parameters
        ---------
        spec_q: Queue-like object
            Spectra are passed to this queue once captured
        capt_q: Queue-like object
            Capture commands are passed to this object using its put() method
        """
        self.continuous_capture = True

        if self.int_time is None:
            raise ValueError('Cannot acquire sequence until initial integration time is correctly set')

        # Setup queue
        capt_q = self._q_check(capt_q, q_type='capt')
        spec_q = self._q_check(spec_q, q_type='spec')

        # Get acquisition rate in seconds
        frame_rep = round(1 / self.framerate)

        # Previous second value for check that we don't take 2 images in one second
        prev_sec = None

        while True:

            # Rethink this later - how to react perhaps depends on what is sent to the queue?
            try:
                mess = capt_q.get(block=False)
                if 'exit_cont' in mess:
                    if mess['exit_cont']:
                        self.continuous_capture = False
                        return

                if 'auto_ss' in mess:
                    # If auto_ss is changed we need to readjust all parameters
                    if not mess['auto_ss']:
                        self.auto_ss = False
                    else:
                        self.auto_ss = True

                    # If we aren't using auto_ss, check for ss in message to set shutter speed
                if not self.auto_ss:
                    if 'ss' in mess:
                        self.int_time = mess['ss']

                if 'framerate' in mess:
                    # We readjust to requested framerate regardless of if auto_ss is True or False
                    frame_rep = round(1 / mess['framerate'])

            except queue.Empty:
                # If there is nothing in the queue telling us to stop then we continue with acquisitions
                pass

            # Get current time
            time_obj = datetime.datetime.now()

            # Only capture an image if we are at the right time
            if time_obj.second % frame_rep == 0 and time_obj != prev_sec:

                # Generate time string
                time_str = format_time(time_obj, self.file_datestr)

                # Acquire spectra
                self.get_spec()

                # Generate filename
                filename = self.generate_filename(time_str, self.file_spec_type['meas'])

                # Add spectrum and filename to queue
                spec_q.put([filename, self.spectrum])

                # Check image saturation and adjust shutter speed if required
                if self.auto_int:
                    adj_saturation = self.check_saturation()
                    if adj_saturation:
                        # Adjust ss_idx, but if we have gone beyond the indices available in ss_list it will throw an
                        # idx error, so we catch this and continue with same int if there are no higher/lower options
                        try:
                            self.int_time_idx += adj_saturation
                            # self.int_time = self.int_list[self.int_time_idx]
                        except IndexError:
                            pass

                # Set seconds value (used as check to prevent 2 images being acquired in same second)
                prev_sec = time_obj.second

    def capture_darks(self):
        """Capture dark images from all shutter speeds in <self.ss_list>"""
        # Loop through shutter speeds in ss_list
        for int_time in self.int_list:

            # Set camera shutter speed
            self.int_time = int_time

            # Get time for stamping
            time_str = format_time(datetime.datetime.now(), self.file_datestr)

            # Acquire image
            self.get_spec()

            # Generate filename for spectrum
            filename = self.generate_filename(time_str, self.file_spec_type['dark'])

            # Add data to queue
            self.spec_q.put(filename)
            self.spec_q.put(self.spectrum)


class SpectrometerConnectionError(Exception):
    """
    Error raised if no spectrometer is detected
    """
    pass


class Scanner:
    """
    Controls arduino to move scanner head
    """
    def __init__(self, ard_com):
        self.ard_com = ard_com                      # COM port for arduino
        self.baud_rate = 9600                       # Baud rate
        self.scan_properties = ScanProperties()     # Object containing important scan properties
        self.arduino = None                         # Arduino serial object

        # Connect to arduino via serial interface
        self.connect_ard()

    def connect_ard(self):
        """Connects to arduino"""
        # Setup arduino serial port
        if self.ard_com is not None:
            try:
                self.arduino = serial.Serial(self.ard_com, self.baud_rate)
            except serial.serialutil.SerialException:
                print('Could not open port to arduino, please check connection and restart program')
        else:
            self.arduino = None
            print('No arduino COM port specified')

    def step_forward(self):
        """Steps the motor forward one step"""
        self.arduino.write(self.scan_properties.scan_fwd)
        reply = self.arduino.read()

    def step_backward(self):
        """Steps the motor forward one step"""
        self.arduino.write(self.scan_properties.scan_back)
        reply = self.arduino.read()

    def scan_return(self):
        """Steps the motor back to clear sky position 1, assuming that the motor has scanned to the end of total_scan_
        steps"""
        for i in range(self.scan_properties.total_scan_steps):
            self.step_backward()

    def full_return(self):
        """Steps scanner all the way back to the start"""
        for i in range(self.scan_properties.step_360):
            self.step_backward()

    def full_forwards(self):
        """Steps scanner all the way to the end"""
        for i in range(self.scan_properties.step_360):
            self.step_forward()

    def clear_to_dark(self):
        """Steps into dark position assuming the the motor is currently at the starting clear-sky position"""
        for i in range(self.scan_properties.dark_steps):
            self.step_backward()

    def dark_to_clear(self):
        """Steps from dark position back to the first clear-sky position"""
        for i in range(self.scan_properties.dark_steps):
            self.step_forward()


class ScanProperties:
    """
    Handles overall control of spectral acquisition, including motor movement
    """
    def __init__(self):
        self.scan_incr = 1.8  # Stepper motors scan increment
        self.scan_range_full = 100  # Range of motion in scanner
        self.total_scan_steps = int(
            np.ceil(self.scan_range_full / self.scan_incr))     # Max num of steps to reset motor at start
        self.step_360 = int(np.ceil(360/self.scan_incr))        # Full number of steps to perform 360 revolution
        self.scan_fwd = b'\x00'  # Set here whether 0 or 1 steps the scanner in the forward direction
        self.scan_back = b'\x01'

        self.dark_steps = 60    # Number of steps to get to dark position
