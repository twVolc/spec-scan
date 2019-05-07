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

class AcquisitionFrame:
    """
    Frame for controlling acquisition settings and instigating acquisitions
    """
    def __init__(self, frame):
        self.setts = SettingsGUI()      # Import settings


        self.start_int_time = 100       # Integration time to load program with
        self.start_scan_range = 45      # Scan angle range to load program with
        try:
            self.spec_ctrl = SpecCtrl(int_time=self.start_int_time)  # Holds spectral control class
        except SpectrometerConnectionError:
            print('Warning!!! No spectrometer detected. Please connect now.')
            self.spec_ctrl = None

        # Setup main frame
        self.frame = ttk.LabelFrame(frame, text='Acquisition Settings', relief=tk.GROOVE, borderwidth=5)

        row = 0     # Row will be incremented for easy griding

        # Single or scan-mode
        self.acq_mode = tk.IntVar()
        self.acq_mode.set(1)
        # self.single_mode = tk.Radiobutton(self.frame, variable=self.acq_mode, text='Single-mode', value=0,
        #                                   bg=self.setts.bgColour,
        #                                   font=self.setts.mainFont).grid(row=row, column=0, sticky='w',
        #                                                                  padx=self.setts.px, pady=self.setts.py)
        # self.scan_mode = tk.Radiobutton(self.frame, variable=self.acq_mode, text='Scan-mode', value=1,
        #                                   bg=self.setts.bgColour,
        #                                 font=self.setts.mainFont).grid(row=row, column=1, sticky='w',
        #                                                                padx=self.setts.px, pady=self.setts.py)
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

        # Acquisition button
        self.acquire_butt = ttk.Button(self.frame, text='ACQUIRE', command=self.acquisition_handler,
                                       ).grid(row=row, column=0, columnspan=2, sticky='nsew')

    def acquisition_handler(self):
        """Controls acquisitions"""
        # If we don't already have a spectrometer set, search for it. If we can't find one, return.
        if self.spec_ctrl is None:
            try:
                self.spec_ctrl = SpecCtrl(int_time=self.start_int_time)  # Holds spectral control class
            except SpectrometerConnectionError:
                print('Warning!!! No spectrometer detected. Please connect now.')
                return

        if self.spec_ctrl.spec is None:
            try:
                self.spec_ctrl.find_device()
            except SpectrometerConnectionError:
                print('No spectrometer found')
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

