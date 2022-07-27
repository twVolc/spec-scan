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
# from doas_routine import DOASWorker, ScanProcess
# from save_spec import SaveSpectra

import numpy as np
from datetime import datetime
import time
import os
import serial
import threading
import queue


class AcquisitionFrame:
    """
    Frame for controlling acquisition settings and instigating acquisitions
    This class brings together the the DOAS work to control processing and plotting of DOAS too
    """
    def __init__(self, frame, doas_worker=None, scan_proc=None,
                 spec_plot=None, doas_plot=None, cd_plot=None, ard_com=None, save_path='C:\\'):
        self.scan_cont = ScanProperties()
        self.scan_proc = scan_proc

        self.setts = SettingsGUI()      # Import settings
        self.doas_worker = doas_worker  # Setup DOASWorker object, used for processing
        self.spec_plot = spec_plot      # Setup SpectraPlot object, used for plotting spectra
        self.doas_plot = doas_plot
        self.cd_plot = cd_plot
        self.save_obj = self.doas_worker.save_obj

        # PROBABLY SETUP THIS PATH THROUGH A FUNCTION WHICH CREATES A NEW DATE DIRECTORY
        # OR THIS MAY BE SETUP OUTSIDE OF THIS CLASS - BY THE MAIN CLASS
        self.save_path_root = save_path     # Root save path (not including the date directory)
        self.save_path = None               # Save path including the date directory
        self.save_path_ind = None           # Save path for individual spectrum acquisitions
        self.scan_dir = None                # Save path for scan
        self.str_len_max = 28
        self.dark_path = None
        self.clear_path = None
        self.plume_path = None
        self.doas_path = None

        self.scanning = False
        self.scan_q = queue.Queue()

        self.start_int_time = 100       # Integration time to load program with
        self.start_coadd = 1
        self.start_scan_range = 45      # Scan angle range to load program with

        # Try to initiate spectrometer control class
        try:
            self.spec_ctrl = SpecCtrl(int_time=self.start_int_time, coadd=self.start_coadd)  # Holds spectral control class

            # Update save object with spec ctrl if we can initiate it
            self.save_obj.spec_ctrl = self.spec_ctrl
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
            print('No arduino COM port specified')

        # Instigates the setup of the GUI
        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):

        # Setup main frame
        self.frame = tk.LabelFrame(frame, text='Acquisition Settings', relief=tk.RAISED, borderwidth=5)

        row = 0 # Row will be incremented for easy griding


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

        # Coadding
        label = ttk.Label(self.frame, text='Co-adding spectra:').grid(row=row, column=0, sticky='w',
                                                                      padx=self.setts.px, pady=self.setts.py)
        self.coadd = tk.IntVar()
        self.coadd.set(self.start_coadd)
        self.coadd_box = tk.Spinbox(self.frame, from_=1, to=100, increment=1, width=3,
                                    textvariable=self.coadd, command=self.update_coadd)
        self.coadd_box.grid(row=row, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)

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

        label = ttk.Label(self.frame2, text='Save directory:').grid(row=0, column=0, sticky='e',
                                                                        padx=self.setts.px, pady=self.setts.py)
        self.save_label = ttk.Label(self.frame2, text=self.save_path_root, width=25)
        self.save_label.grid(row=0, column=1, sticky='nsew', padx=self.setts.px, pady=self.setts.py)
        self.__config_save_path__(self.save_path_root)
        self.save_butt = ttk.Button(self.frame2, text='Choose directory', command=self.change_save_path)
        self.save_butt.grid(row=1, column=1, sticky='nsew', padx=self.setts.px, pady=self.setts.py)

        self.label3 = ttk.Label(self.frame2, text='Dark spectrum file:').grid(row=2, column=0, sticky='w',
                                                                             padx=self.setts.px, pady=self.setts.py)
        self.dark_file_label = ttk.Label(self.frame2, text='N/A', width=25)
        self.dark_file_label.grid(row=2, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.dark_button = ttk.Button(self.frame2, text='Dark Capture',
                                      command=self.dark_capture).grid(row=3, column=1, sticky='nsew')
        self.label4 = ttk.Label(self.frame2, text='Clear spectrum file:').grid(row=4, column=0, sticky='w',
                                                                              padx=self.setts.px, pady=self.setts.py)
        self.clear_file_label = ttk.Label(self.frame2, text='N/A', width=25)
        self.clear_file_label.grid(row=4, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.clear_button = ttk.Button(self.frame2, text='Clear Capture',
                                       command=self.clear_capture).grid(row=5, column=1, sticky='nsew')

        row += 1

        # Acquisition button
        self.acquire_butt = ttk.Button(self.frame, text='ACQUIRE', command=self.acquisition_handler)
        self.acquire_butt.grid(row=row, column=0, columnspan=2, sticky='nsew')

        self.stop_acquire_butt = ttk.Button(self.frame, text='STOP SCAN!', command=self.stop_scan)
        self.stop_acquire_butt.grid(row=row, column=0, columnspan=2, sticky='nsew')
        self.stop_acquire_butt.grid_remove()

    def change_save_path(self):
        """Brings up dialog box to change the save directory"""
        save_path = filedialog.askdirectory(initialdir=self.save_path_root, title='Select save directory')
        if save_path:
            save_path += '/'
            self.__config_save_path__(save_path)

    def __config_save_path__(self, save_path):
        """Configures label and save paths"""
        self.save_path = save_path
        if len(save_path) > self.str_len_max:
            save_path = '...' + self.save_path_root[-(self.str_len_max-5):]
        self.save_label.configure(text=save_path)

        # Make date directory to save spectra in
        self.save_path = self.save_path_root + datetime.now().strftime('%Y-%m-%d') + '/'
        if not os.path.exists(self.save_path):
            os.mkdir(self.save_path)

        # Make save path for individual acquisitions
        self.save_path_ind = self.save_path + 'Test_spectra/'
        if not os.path.exists(self.save_path_ind):
            os.mkdir(self.save_path_ind)

    def update_coadd(self):
        """Update coadding value"""
        try:
            self.spec_ctrl.coadd = self.coadd.get()
        # If no spectrometer is detected we just ignore the error thrown
        except AttributeError:
            pass

    def dark_capture(self):
        """Controls dark spectrum capture"""
        if not self._check_connection():
            return

        # Set integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Set doas_worker wavelength attribute to the wavelength calibration of spectrometer
        self.doas_worker.wavelengths = self.spec_ctrl.wavelengths

        # Step to dark position
        for i in range(self.scan_cont.dark_steps):
            # Step motor
            self.arduino.write(self.scan_cont.scan_back)
            reply = self.arduino.read()

        # Ignore first spectrum as it could have been acquired prior to
        self.doas_worker.dark_spec = self.spec_ctrl.get_spec()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        dark_filename = '{}_dark.txt'.format(time)
        if self.scanning:
            self.dark_path = self.scan_dir + dark_filename
        elif not self.scanning:
            self.dark_path = self.save_path_ind + dark_filename
        else:
            print('Unrecognised scanning flag. Stopping acquisition')
            return
        # self.doas_worker.save_dark(self.dark_path)
        self.save_obj.save_dark(self.dark_path)

        # Change GUI label for dark file
        self.dark_file_label.configure(text=dark_filename)

        # Update plot with new data
        self.spec_plot.update_dark()

        # Step back to starting position
        for i in range(self.scan_cont.dark_steps):
            # Step motor
            self.arduino.write(self.scan_cont.scan_fwd)
            reply = self.arduino.read()

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
        if self.scanning:
            self.clear_path = self.scan_dir + clear_filename
        elif not self.scanning:
            self.clear_path = self.save_path_ind + clear_filename
        else:
            print('Unrecognised scanning flag. Stopping acquisition')
            return
        # self.doas_worker.save_clear_raw(self.clear_path)
        self.save_obj.save_clear_raw(self.clear_path)

        # Change GUI label for dark file
        self.clear_file_label.configure(text=clear_filename)

        # Update plot with new data
        self.spec_plot.update_clear()

    def save_processed_spec(self):
        """Saves processed spectrum with all useful information"""
        if self.doas_worker.wavelengths_cut is None or self.doas_worker.abs_spec_cut is None:
            print('No processing information, save not completed')
            return
        if self.doas_worker.ref_spec_types[0] not in self.doas_worker.ref_spec_fit:
            print('No processing information, save not completed')
            return

        # Setup directory and filename
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        doas_filename = '{}_doas.txt'.format(time)
        if self.scanning:
            doas_dir = self.scan_dir + 'Processing_1/'
        else:
            doas_dir = self.save_path_ind + 'Processing_1/'
        if not os.path.exists(doas_dir):
            os.mkdir(doas_dir)
        self.doas_path = doas_dir + doas_filename

        print(self.doas_path)

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
                          self.doas_worker.column_density['SO2']))

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
        if self.scanning:
            self.plume_path = self.scan_dir + plume_filename
        elif not self.scanning:
            self.plume_path = self.save_path_ind + plume_filename
        else:
            print('Unrecognised scanning flag. Stopping acquisition')
            return
        # self.doas_worker.save_plume_raw(self.plume_path)
        self.save_obj.save_plume_raw(self.plume_path)

        # Update plot with new data
        self.spec_plot.update_plume()

        # Try processing data
        self.doas_worker.process_doas()

        # Update plot and save processed spectrum if data was processed
        if self.doas_worker.processed_data:
            self.doas_plot.update_plot()
            # self.save_processed_spec()

            # Setup directory and filename
            time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
            doas_filename = '{}_doas.txt'.format(time)
            if self.scanning:
                doas_dir = self.scan_dir + 'Processing_1/'
            else:
                doas_dir = self.save_path_ind + 'Processing_1/'
            if not os.path.exists(doas_dir):
                os.mkdir(doas_dir)
            self.doas_path = doas_dir + doas_filename

            self.save_obj.save_processed_spec(self.dark_path, self.clear_path, self.plume_path, self.doas_path)

    def acquisition_handler(self):
        """Controls acquisitions"""
        if not self._check_connection():
            return

        # Implement function determined by radiobutton decision
        acq_mode = self.acq_mode.get()
        if acq_mode == 0:
            self.plume_capture()
        elif acq_mode == 1:
            self.acquire_butt.grid_remove()
            self.stop_acquire_butt.grid()
            self.scan_thread = threading.Thread(target=self.acquire_scan)
            self.scan_thread.start()
            # self.acquire_scan()
        else:
            raise ValueError('Acquisition mode expected 0 or 1. Got {}'.format(acq_mode))

    def stop_scan(self):
        """Stops scanning"""
        # Instigate scan stop
        self.scan_q.put(1)

        # # Wait for scan thread to finish
        # self.scan_thread.join()

        # Update buttons to allow scanning requests again
        self.stop_acquire_butt.grid_remove()
        self.acquire_butt.grid()

    def acquire_scan(self):
        """Perform scan acquisition"""
        self.scanning = True

        # Update plume distance and speed from GUI
        self.__update_scan_params__()

        # Calculate number of steps needed to perform requested scan
        num_steps = int(self.scan_range.get() / self.scan_cont.scan_incr)

        while True:
            # Setup scan directory
            self.scan_dir = self.save_path + 'Scan_1/'
            idx = 2
            while os.path.exists(self.scan_dir):
                self.scan_dir = self.save_path + 'Scan_{}/'.format(idx)
                idx += 1
            os.mkdir(self.scan_dir)

            # Clear any previous plots of emissions rates
            self.cd_plot.clear_plot()

            # Clear scan data so that fresh arrays are present
            self.scan_proc.clear_data()

            # Acquire dark image
            self.dark_capture()

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
                self.scan_proc.add_data(scan_angle, self.doas_worker.column_density['SO2'])

                # Update column densities plot
                self.cd_plot.update_plot()

                # # Sleep to attempt to give GUI time to update everything
                # time.sleep(0.25)

            self.scan_proc.calc_emission_rate()
            self.cd_plot.update_emission_rate()
            # self.save_scan()
            self.save_obj.save_scan(self.scan_dir + 'Processing_1/')

            # Reset stepper motor
            for i in range(num_steps):
                self.arduino.write(self.scan_cont.scan_back)
                reply = self.arduino.read()

            # Check to see if we have been sent a message to stop scanning
            try:
                quit = self.scan_q.get(block=False)
                # Clear queue in case we sent multiple commands to stop
                with self.scan_q.mutex:
                    self.scan_q.queue.clear()
                self.scanning = False
                return
            except queue.Empty:
                continue

    def save_scan(self):
        """Saves scan information"""
        filename = 'Scan_data.txt'
        doas_dir = self.scan_dir + 'Processing_1/'
        try:
            np.savetxt(doas_dir + filename, np.transpose([self.scan_proc.scan_angles,
                                                          self.scan_proc.column_densities]),
                        header='Processed scan details\n'
                               'Plume speed: {}\n'
                               'Plume distance: {}\n'
                               'Emission rate [kg/s]: {}\n'
                               'Emission rate [t/day]: {}\n'
                               'Scan angle [deg]\tColumn density [ppm.m]'.format(
                            self.scan_proc.plume_speed, self.scan_proc.plume_distance,
                            self.scan_proc.SO2_flux, self.scan_proc.flux_tons))
        except Exception as e:
            print(e)

    def __update_scan_params__(self):
        """Updates scan paramters"""
        self.scan_proc.plume_speed = self.cd_plot.plume_speed
        self.scan_proc.plume_distance = self.cd_plot.plume_dist

    def _check_connection(self):
        """Checks spectrometer connection"""
        # If we don't already have a spectrometer set, search for it. If we can't find one, return.
        if self.spec_ctrl is None:
            try:
                self.spec_ctrl = SpecCtrl(int_time=self.start_int_time, coadd=self.coadd.get())  # Holds spectral control class

                # Update save object with spec ctrl if we can initiate it
                self.save_obj.spec_ctrl = self.spec_ctrl
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
