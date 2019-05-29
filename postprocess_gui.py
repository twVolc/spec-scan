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
from plotting_gui import SpectraPlot
from doas_routine import DOASWorker, ScanProcess

import numpy as np
from datetime import datetime

class PostProcess:
    """
    Creates frame for postprocessing and handles procedures
    """
    def __init__(self, frame, doas_worker=DOASWorker(), scan_proc=ScanProcess(),
                 spec_plot=None, doas_plot=None, cd_plot=None, init_dir='C:\\'):
        self.doas_worker = doas_worker
        self.scan_proc = scan_proc
        self.spec_plot = spec_plot
        self.doas_plot = doas_plot
        self.cd_plot = cd_plot

        self.str_len_max = 30
        self.pdx = 5
        self.pdy = 5

        self.init_dir = init_dir    # Directory to begin filedialog prompts
        self.save_path = init_dir   # Directory to save processed date to
        self.dark_path = None
        self.clear_path = None
        self.plume_path = None
        self.doas_path = None

        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):
        self.frame = ttk.LabelFrame(frame, text='Post-Processing', relief=tk.RAISED, borderwidth=5)

        # Dark load
        label = ttk.Label(self.frame, text='Dark spectrum file:')
        self.dark_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_dark_butt = ttk.Button(self.frame, text='Load Dark', command=self.load_dark)

        row = 0
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.dark_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_dark_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)

        # Clear load
        label = ttk.Label(self.frame, text='Clear spectrum file:')
        self.clear_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_clear_butt = ttk.Button(self.frame, text='Load Clear', command=self.load_clear)

        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.clear_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_clear_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)

        # Plume load
        label = ttk.Label(self.frame, text='Dark spectrum file:')
        self.plume_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_plume_butt = ttk.Button(self.frame, text='Load Plume', command=self.load_plume)

        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.plume_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_plume_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)

        # Load scan
        label = ttk.Label(self.frame, text='Scan directory:')
        self.scan_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_scan_butt = ttk.Button(self.frame, text='Load Scan', command=self.load_scan)

        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.scan_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_scan_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)


    def load_dark(self):
        """Load dark spectrum"""
        # Bring up dialog to find file
        self.dark_path = filedialog.askopenfilename(initialdir=self.init_dir,
                                                    title='Select dark spectrum file',
                                                    filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        if not self.dark_path:
            return

        self.init_dir, dark_filename = self.dark_path.rsplit('/', maxsplit=1)

        # Update label in widget
        if len(dark_filename) > self.str_len_max:
            self.dark_label.configure(text='...' + dark_filename[-(self.str_len_max - 3):])
        else:
            self.dark_label.configure(text=dark_filename)

        # Extract data
        data = np.loadtxt(self.dark_path)
        self.doas_worker.wavelengths, self.doas_worker.dark_spec = data.T

        # Update dark plot
        self.spec_plot.update_dark()

    def load_clear(self):
        """Load clear spectrum"""
        # Bring up dialog to find file
        self.clear_path = filedialog.askopenfilename(initialdir=self.init_dir,
                                                    title='Select clear spectrum file',
                                                    filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        if not self.clear_path:
            return

        self.init_dir, clear_filename = self.clear_path.rsplit('/', maxsplit=1)
        # Update label in widget
        if len(clear_filename) > self.str_len_max:
            self.clear_label.configure(text='...' + clear_filename[-(self.str_len_max - 3):])
        else:
            self.clear_label.configure(text=clear_filename)

        # Extract data
        data = np.loadtxt(self.clear_path)
        self.doas_worker.wavelengths, self.doas_worker.clear_spec_raw = data.T

        # Update dark plot
        self.spec_plot.update_clear()

    def load_plume(self):
        """Load plume spectrum"""
        # Bring up dialog to find file
        self.plume_path = filedialog.askopenfilename(initialdir=self.init_dir,
                                                    title='Select in-plume spectrum file',
                                                    filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        if not self.plume_path:
            return

        self.init_dir, plume_filename = self.plume_path.rsplit('/', maxsplit=1)
        # Update label in widget
        if len(plume_filename) > self.str_len_max:
            self.plume_label.configure(text='...' + plume_filename[-(self.str_len_max - 3):])
        else:
            self.plume_label.configure(text=plume_filename)

        # Extract data
        data = np.loadtxt(self.plume_path)
        self.doas_worker.wavelengths, self.doas_worker.plume_spec_raw = data.T

        # Update dark plot
        self.spec_plot.update_plume()

        # Try processing data
        self.doas_worker.process_doas()

        # Update doas plot
        self.doas_plot.update_plot()

        # Save processed spectrum if data was processed
        if self.doas_worker.processed_data:
            self.save_processed_spec()

    def load_scan(self):
        """Load scan directory"""
        pass

    def save_processed_spec(self):
        """Saves processed spectrum with all useful information"""
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        doas_filename = '{}_doas.txt'.format(time)
        self.doas_path = self.save_path + doas_filename

        np.savetxt(self.doas_path, np.transpose([self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit['SO2'],
                                                 self.doas_worker.abs_spec_cut]),
                   header='Processed DOAS spectrum\n'
                          'Dark spectrum: {}\nClear spectrum{}\nPlume spectrum{}\n'
                          'Shift: {}\nStretch: {}\n'
                          'Stray range [nm]: {}:{}\nFit window [nm]: {}:{}\n'
                          'Column density [ppm.m]: {}\n'
                          'Wavelength [nm]\tReference spectrum (fitted)\tAbsorbance spectrum'.format(
                          self.dark_path, self.clear_path, self.plume_path,
                          self.doas_worker.shift, self.doas_worker.stretch,
                          self.doas_worker.start_stray_wave, self.doas_worker.end_stray_wave,
                          self.doas_worker.start_fit_wave, self.doas_worker.end_fit_wave,
                          self.doas_worker.column_amount))