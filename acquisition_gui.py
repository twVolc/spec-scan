import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from tkinter import messagebox
import tkinter.font as tkFont

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from gui_subs import SettingsGUI
from controllers import SpecCtrl, SpectrometerConnectionError, ScanProperties
from plotting_gui import SpectraPlot
from doas_routine import DOASWorker, ScanProcess

import numpy as np
from datetime import datetime
import time
import serial


class AcquisitionFrame:
    """
    Frame for controlling acquisition settings and instigating acquisitions
    This class brings together the the DOAS work to control processing and plotting of DOAS too
    """
    def __init__(self, frame, doas_worker=DOASWorker(), scan_proc=ScanProcess(),
                 spec_plot=None, doas_plot=None, cd_plot=None, ard_com=None):
        self.scan_cont = ScanProperties()
        self.scan_proc = scan_proc

        self.setts = SettingsGUI()      # Import settings
        self.doas_worker = doas_worker  # Setup DOASWorker object, used for processing
        self.spec_plot = spec_plot      # Setup SpectraPlot object, used for plotting spectra
        self.doas_plot = doas_plot
        self.cd_plot = cd_plot

        # PROBABLY SETUP THIS PATH THROUGH A FUNCTION WHICH CREATES A NEW DATE DIRECTORY
        # OR THIS MAY BE SETUP OUTSIDE OF THIS CLASS - BY THE MAIN CLASS
        self.save_path = 'C:\\Users\\tw9616\\Documents\\PostDoc\\Scanning Spectrometer\\SpecScan\\Spectra\\'
        self.dark_path = None
        self.clear_path = None
        self.plume_path = None
        self.doas_path = None

        self.start_int_time = 100       # Integration time to load program with
        self.start_scan_range = 45      # Scan angle range to load program with
        try:
            self.spec_ctrl = SpecCtrl(int_time=self.start_int_time)  # Holds spectral control class
        except SpectrometerConnectionError:
            print('Warning!!! No spectrometer detected. Please connect now.')
            self.spec_ctrl = None

        # Setup arduino serial port
        if ard_com is not None:
            try:
                self.arduino = serial.Serial(ard_com, 9600)
            except serial.serialutil.SerialException:
                print('Could not open port to arduino, please check connection and restart program')
        else:
            self.arduino = None

        # Instigates the setup of the GUI
        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):

        # Setup main frame
        self.frame = ttk.LabelFrame(frame, text='Acquisition Settings', relief=tk.RAISED, borderwidth=5)

        row = 0     # Row will be incremented for easy griding

        # Single or scan-mode
        self.acq_mode = tk.IntVar()
        self.acq_mode.set(1)
        self.single_mode = ttk.Radiobutton(self.frame, variable=self.acq_mode, text='Single-mode', value=0,
                                           ).grid(row=row, column=0, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.scan_mode = ttk.Radiobutton(self.frame, variable=self.acq_mode, text='Scan-mode', value=1,
                                         ).grid(row=row, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)

        row += 1

        # Integration time setup
        self.label1 = ttk.Label(self.frame, text='Integration time (ms):').grid(row=row, column=0, sticky='w',
                                                                                padx=self.setts.px, pady=self.setts.py)
        self.int_time = tk.DoubleVar()     # Integration time holder
        self.int_time.set(self.start_int_time)
        self.int_entry = ttk.Entry(self.frame, textvariable=self.int_time, width=5,
                                   ).grid(row=row, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        row += 1

        # Scan options
        self.label2 = ttk.Label(self.frame, text='Scanning range (°):').grid(row=row, column=0, sticky='w',
                                                                             padx=self.setts.px, pady=self.setts.py)
        self.scan_range = tk.DoubleVar()  # Integration time holder
        self.scan_range.set(self.start_scan_range)
        self.scan_entry = ttk.Entry(self.frame, textvariable=self.scan_range, width=5,
                                    ).grid(row=row, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        row += 1

        # Dark and clear (fraunhofer) spectrum acquisitions
        self.frame2 = ttk.Frame(self.frame, relief=tk.GROOVE, borderwidth=5)
        self.frame2.grid(row=row, column=0, columnspan=2, sticky='nsew')
        self.label3 = ttk.Label(self.frame2, text='Dark spectrum file:').grid(row=0, column=0, sticky='w',
                                                                             padx=self.setts.px, pady=self.setts.py)
        self.dark_file_label = ttk.Label(self.frame2, text='N/A', width=25)
        self.dark_file_label.grid(row=0, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.dark_button = ttk.Button(self.frame2, text='Dark Capture',
                                      command=self.dark_capture).grid(row=1, column=1, sticky='nsew')
        self.label4 = ttk.Label(self.frame2, text='Clear spectrum file:').grid(row=2, column=0, sticky='w',
                                                                              padx=self.setts.px, pady=self.setts.py)
        self.clear_file_label = ttk.Label(self.frame2, text='N/A', width=25)
        self.clear_file_label.grid(row=2, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.clear_button = ttk.Button(self.frame2, text='Clear Capture',
                                       command=self.clear_capture).grid(row=3, column=1, sticky='nsew')

        row += 1

        # Acquisition button
        self.acquire_butt = ttk.Button(self.frame, text='ACQUIRE', command=self.acquisition_handler,
                                       ).grid(row=row, column=0, columnspan=2, sticky='nsew')


    def dark_capture(self):
        """Controls dark spectrum capture"""
        if not self._check_connection():
            return

        # Set integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Set doas_worker wavelength attribute to the wavelength calibration of spectrometer
        self.doas_worker.wavelengths = self.spec_ctrl.wavelengths

        # Ignore first spectrum as it could have been acquired prior to
        self.doas_worker.dark_spec = self.spec_ctrl.get_spec()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        dark_filename = '{}_dark.txt'.format(time)
        self.dark_path = self.save_path + dark_filename
        self.doas_worker.save_dark(self.dark_path)

        # Change GUI label for dark file
        self.dark_file_label.configure(text=dark_filename)

        # Update plot with new data
        self.spec_plot.update_dark()


    def clear_capture(self):
        """Controls clear spectrum capture"""
        if not self._check_connection():
            return

        # Set integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Set doas_worker wavelength attribute to the wavelength calibration of spectrometer
        self.doas_worker.wavelengths = self.spec_ctrl.wavelengths

        # Ignore first spectrum as it could have been acquired prior to
        self.doas_worker.clear_spec_raw = self.spec_ctrl.get_spec()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        clear_filename = '{}_clear.txt'.format(time)
        self.clear_path = self.save_path + clear_filename
        self.doas_worker.save_clear_raw(self.clear_path)

        # Change GUI label for dark file
        self.clear_file_label.configure(text=clear_filename)

        # Update plot with new data
        self.spec_plot.update_clear()

    def save_processed_spec(self):
        """Saves processed spectrum with all useful information"""
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        doas_filename = '{}_doas.txt'.format(time)
        self.doas_path = self.save_path + doas_filename

        np.savetxt(self.doas_path, np.transpose([self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit['SO2'],
                                                 self.doas_worker.abs_spec_cut]),
                   header='Processed DOAS spectrum\n'
                          'Dark spectrum: {}\nClear spectrum: {}\nPlume spectrum: {}\n'
                          'Shift: {}\nStretch: {}\n'
                          'Stray range [nm]: {}:{}\nFit window [nm]: {}:{}\n'
                          'Column density [ppm.m]: {}\n'
                          'Wavelength [nm]\tReference spectrum (fitted)\tAbsorbance spectrum'.format(
                          self.dark_path, self.clear_path, self.plume_path,
                          self.doas_worker.shift, self.doas_worker.stretch,
                          self.doas_worker.start_stray_wave, self.doas_worker.end_stray_wave,
                          self.doas_worker.start_fit_wave, self.doas_worker.end_fit_wave,
                          self.doas_worker.column_amount))

    def acquisition_handler(self):
        """Controls acquisitions"""
        if not self._check_connection():
            return

        # Implement function determined by radiobutton decision
        acq_mode = self.acq_mode.get()
        if acq_mode == 0:
            self.plume_capture()
        elif acq_mode == 1:
            self.acquire_scan()
        else:
            raise ValueError('Acquisition mode expected 0 or 1. Got {}'.format(acq_mode))

    def plume_capture(self):
        """Perform single acquisition"""
        # Could possibly condense this work into a a function which is called by plume_capture and acquire scan, as
        # much of this function is used in acquire scan

        # Set spectrometer control object to correct integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Set doas_worker wavelength attribute to the wavelength calibration of spectrometer
        self.doas_worker.wavelengths = self.spec_ctrl.wavelengths

        # Acquire spectrum (this method includes discarding the first spectrum retrieved from the spectrometer)
        self.doas_worker.plume_spec_raw = self.spec_ctrl.get_spec()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        plume_filename = '{}_plume.txt'.format(time)
        self.plume_path = self.save_path + plume_filename
        self.doas_worker.save_plume_raw(self.plume_path)

        # Update plot with new data
        self.spec_plot.update_plume()

        # Try processing data
        self.doas_worker.process_doas()

        # Update plot and save processed spectrum if data was processed
        if self.doas_worker.processed_data:
            self.doas_plot.update_plot()
            self.save_processed_spec()

    def acquire_scan(self):
        """Perform scan acquisition"""
        num_steps = int(self.scan_range.get() / self.scan_cont.scan_incr)

        self.scan_proc.scan_angles = np.array([])
        self.scan_proc.column_densities = np.array([])

        for i in range(self.scan_cont.return_steps):
            self.arduino.write(self.scan_cont.scan_back)
            reply = self.arduino.read()

        # Set clear spectrum as first spectrum acquired
        self.clear_capture()

        scan_angle = 0
        for i in range(num_steps):
            # Step motor
            self.arduino.write(self.scan_cont.scan_fwd)
            reply = self.arduino.read()
            scan_angle += self.scan_cont.scan_incr

            # Wait to make sure stepper has moved
            time.sleep(0.5)

            # Acquire spectrum - includes plot updates and processing
            self.plume_capture()

            # Update scan processing object with new data
            self.scan_proc.add_data(scan_angle, self.doas_worker.column_amount)

            # Update column densities plot
            self.cd_plot.update_plot()

        # Reset stepper motor
        for i in range(num_steps):
            self.arduino.write(self.scan_cont.scan_back)
            reply = self.arduino.read()


    def _check_connection(self):
        """Checks spectrometer connection"""
        # If we don't already have a spectrometer set, search for it. If we can't find one, return.
        if self.spec_ctrl is None:
            try:
                self.spec_ctrl = SpecCtrl(int_time=self.start_int_time)  # Holds spectral control class
            except SpectrometerConnectionError:
                print('Warning!!! No spectrometer detected. Please connect now.')
                return 0

        if self.spec_ctrl.spec is None:
            try:
                self.spec_ctrl.find_device()
            except SpectrometerConnectionError:
                print('No spectrometer found')
                return 0

        return 1    # Returns if connections were successful
