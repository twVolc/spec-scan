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
from controllers import SpecCtrl, SpectrometerConnectionError

from datetime import datetime

class AcquisitionFrame:
    """
    Frame for controlling acquisition settings and instigating acquisitions
    """
    def __init__(self, frame, doas_worker):
        self.setts = SettingsGUI()      # Import settings
        self.doas_worker = doas_worker  # Setup DOASWorker object, used for processing

        # PROBABLY SETUP THIS PATH THROUGH A FUNCTION WHICH CREATES A NEW DATE DIRECTORY
        # OR THIS MAY BE SETUP OUTSIDE OF THIS CLASS - BY THE MAIN CLASS
        self.save_path = 'C:\\Users\\tw9616\\Documents\\PostDoc\\Scanning Spectrometer\\SpecScan\\Spectra\\'
        self.dark_filename = None
        self.clear_filename = None
        self.plume_filename = None

        self.start_int_time = 100       # Integration time to load program with
        self.start_scan_range = 45      # Scan angle range to load program with
        try:
            self.spec_ctrl = SpecCtrl(int_time=self.start_int_time)  # Holds spectral control class
        except SpectrometerConnectionError:
            print('Warning!!! No spectrometer detected. Please connect now.')
            self.spec_ctrl = None

        # Instigates the setup of the GUI
        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):

        # Setup main frame
        self.frame = ttk.LabelFrame(frame, text='Acquisition Settings', relief=tk.GROOVE, borderwidth=5)

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
        self.dark_file_label = ttk.Label(self.frame2, text='N/A')
        self.dark_file_label.grid(row=0, column=1, sticky='w', padx=self.setts.px, pady=self.setts.py)
        self.dark_button = ttk.Button(self.frame2, text='Dark Capture',
                                      command=self.dark_capture).grid(row=1, column=1, sticky='nsew')
        self.label4 = ttk.Label(self.frame2, text='Clear spectrum file:').grid(row=2, column=0, sticky='w',
                                                                              padx=self.setts.px, pady=self.setts.py)
        self.clear_file_label = ttk.Label(self.frame2, text='N/A')
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

        # Set doas_worker wavelength attribute to the wavelength calibration of spectrometer
        self.doas_worker.wavelengths = self.spec_ctrl.wavelengths

        # Ignore first spectrum as it could have been acquired prior to
        self.doas_worker.dark_spec = self.spec_ctrl.get_spec()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        self.dark_filename = '{}_dark.txt'.format(time)
        self.doas_worker.save_dark(self.save_path + self.dark_filename)

        # Change GUI label for dark file
        self.dark_file_label.configure(text=self.dark_filename)


    def clear_capture(self):
        """Controls clear spectrum capture"""
        if not self._check_connection():
            return



    def acquisition_handler(self):
        """Controls acquisitions"""
        if not self._check_connection():
            return

        # Implement function determined by radiobutton decision
        acq_mode = self.acq_mode.get()
        if acq_mode == 0:
            self.acquire_single()
        elif acq_mode == 1:
            self.acquire_scan()
        else:
            raise ValueError('Acquisition mode expected 0 or 1. Got {}'.format(acq_mode))


    def acquire_single(self):
        """Perform single acquisition"""
        # Set spectrometer control object to correct integration time
        self.spec_ctrl.int_time = self.int_time.get()
        spectrum = self.spec_ctrl.get_spec()
        print(spectrum)

    def acquire_scan(self):
        """Perform scan acquisition"""
        # Might want to discard first sepctrum captured at each step, because the spectrometer is constantly acquiring,
        # so the first spectrum may have been acquired whilst the spectrometer was facing a different direction
        pass



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
