import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from datetime import datetime
import numpy as np

from controllers import SpecCtrl, SpectrometerConnectionError
from doas_routine import DOASWorker
from gui_subs import SettingsGUI


class CalPlot:
    """
    Generates plot for calibrating the spectrometer using a mercury lamp
    Instrument line shape is extracted from an emission line of the user's choice
    """
    def __init__(self, frame, spec_ctrl, doas_worker=DOASWorker(), ILS=None):
        self.spec_ctrl = spec_ctrl      # Must be instance of SpecCtrl
        self.doas_worker = doas_worker

        self.save_path = self.save_path = 'C:\\Users\\tw9616\\Documents\\PostDoc\\Scanning Spectrometer\\SpecScan\\Spectra\\Cal\\'
        self.str_len_max = 25   # Maximum string length allowed for filenames before they are abbreviated

        self.start_int_time = 100   # Starting integration time
        self.max_DN = 2**16 - 1

        self.wavelengths = None
        self.dark_spec = None
        self.cal_spec_raw = None
        self.cal_spec_corr = None   # Dark corrected calibration spectrum

        # Setup widget
        self.__setup_gui__(frame)

        # If ILS path is provided for pre-loading, load it
        self.ILS_path = ILS
        if self.ILS_path is not None:
            self.load_ILS()
        else:
            self.ILS_wavelengths = None
            self.ILS = None

    def __setup_gui__(self,frame):

        self.frame = ttk.Frame(frame)

        self.frame_cal = ttk.Frame(self.frame, relief=tk.RAISED, borderwidth=4)
        self.frame_cal.pack(side=tk.TOP)
        self.frame_ILS = ttk.Frame(self.frame, relief=tk.RAISED, borderwidth=4)
        self.frame_ILS.pack(side=tk.RIGHT, anchor='n')

        # -----------------------------------------------
        # WIDGET SETUP
        # ------------------------------------------------
        self.frame2 = ttk.Frame(self.frame_cal)
        self.frame2.pack(side=tk.TOP, fill=tk.X, expand=1)

        # Dark and calibration spectrum acquisitions
        self.frame_acq = ttk.Frame(self.frame2, relief=tk.GROOVE, borderwidth=5)
        self.frame_acq.pack(side=tk.TOP, anchor='w')

        # Acquisition settings
        label = ttk.Label(self.frame_acq, text='Integration time (ms):').grid(row=0, column=0, sticky='w')
        self.int_time = tk.DoubleVar()  # Integration time holder
        self.int_time.set(self.start_int_time)
        self.int_entry = ttk.Entry(self.frame_acq, textvariable=self.int_time, width=5).grid(row=0, column=1, sticky='w')

        label = ttk.Label(self.frame_acq, text='Dark spectrum file:').grid(row=1, column=0, sticky='w')
        self.dark_file_label = ttk.Label(self.frame_acq, text='N/A', width=25)
        self.dark_file_label.grid(row=1, column=1, sticky='w')
        self.dark_button = ttk.Button(self.frame_acq, text='Dark Capture',
                                      command=self.dark_capture).grid(row=1, column=2, sticky='nsew')
        label = ttk.Label(self.frame_acq, text='Calibration spectrum file:').grid(row=2, column=0, sticky='w')
        self.cal_file_label = ttk.Label(self.frame_acq, text='N/A', width=25)
        self.cal_file_label.grid(row=2, column=1, sticky='w')
        self.cal_button = ttk.Button(self.frame_acq, text='Calibration Capture',
                                       command=self.cal_capture).grid(row=2, column=2, sticky='nsew')


        # ILS extraction
        self.frame_lines = ttk.Frame(self.frame_cal, relief=tk.GROOVE, borderwidth=5)
        self.frame_lines.pack(side=tk.TOP, anchor='w')

        self.ILS_start = tk.DoubleVar()
        self.ILS_start.set(self.doas_worker.start_fit_wave)
        self.ILS_box_start = tk.Spinbox(self.frame_lines, from_=0, to=400, increment=0.1, width=5,
                                          textvariable=self.ILS_start, command=self.update_ILS_start)
        self.ILS_end = tk.DoubleVar()
        self.ILS_end.set(self.doas_worker.end_fit_wave)
        self.ILS_box_end = tk.Spinbox(self.frame_lines, from_=1, to=400, increment=0.1, width=5,
                                        textvariable=self.ILS_end, command=self.update_ILS_end)

        label = tk.Label(self.frame_lines, text='Emission line start:').pack(side=tk.LEFT)
        self.ILS_box_start.pack(side=tk.LEFT)
        label = tk.Label(self.frame_lines, text='Emission line end:').pack(side=tk.LEFT)
        self.ILS_box_end.pack(side=tk.LEFT)

        # ------------------------------------------------
        # FIGURE SETUP
        # ------------------------------------------------
        self.fig = plt.Figure(figsize=(10, 3), dpi=100)

        self.ax = self.fig.subplots(1, 1)
        self.ax.set_ylabel('DN')
        self.ax.set_ylim([0, self.max_DN])
        self.ax.set_xlim([250, 400])
        self.ax.set_xlabel('Wavelength [nm]')
        self.ax.grid(True)
        self.plt_colours = ['b', 'r']
        for i in range(2):
            self.ax.plot([250, 400], [0, 0], self.plt_colours[i], linewidth=1)
        self.ax.legend(('Dark', 'Calibration'), loc=1, framealpha=1)
        self.fig.tight_layout()

        # Stray light
        self.min_ILS_line = self.ax.plot([self.doas_worker.start_fit_wave, self.doas_worker.start_fit_wave],
                                           [0, self.max_DN], 'y')
        self.max_ILS_line = self.ax.plot([self.doas_worker.end_fit_wave, self.doas_worker.end_fit_wave],
                                           [0, self.max_DN], 'y')

        width = self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave
        self.ILS_range = Rectangle((self.doas_worker.start_fit_wave, 0), width=width, height=self.max_DN,
                                     fill=True, facecolor='y', alpha=0.2)
        self.ILS_range_patch = self.ax.add_patch(self.ILS_range)

        # Organise and draw canvas
        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame_cal)
        self.canv.draw()
        self.canv.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # --------------------------------------------------------------------------------------------------------------
        # ILS FUCTIONS AND PLOT
        # --------------------------------------------------------------------------------------------------------------
        # FUNCTIONS
        self.frame_ILS_func = ttk.Frame(self.frame_ILS, relief=tk.GROOVE, borderwidth=5)
        self.frame_ILS_func.pack(side=tk.LEFT, anchor='n')

        label1 = tk.Label(self.frame_ILS_func, text='Loaded ILS file:')
        self.ILS_load_label = tk.Label(self.frame_ILS_func, text='N/A', width=self.str_len_max)
        self.load_ILS_button = ttk.Button(self.frame_ILS_func, text='Load ILS', command=self.choose_ILS)

        label2 = tk.Label(self.frame_ILS_func, text='Saved ILS file:')
        self.ILS_save_label = tk.Label(self.frame_ILS_func, text='N/A', width=self.str_len_max)
        self.save_ILS_button = ttk.Button(self.frame_ILS_func, text='Save ILS', command=self.save_ILS)

        label1.grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.ILS_load_label.grid(row=0, column=1, padx=5, pady=5)
        self.load_ILS_button.grid(row=0, column=2, padx=5, pady=5)

        label2.grid(row=1, column=0, sticky='w', padx=5, pady=5)
        self.ILS_save_label.grid(row=1, column=1, padx=5, pady=5)
        self.save_ILS_button.grid(row=1, column=2, padx=5, pady=5)




        # PLOT
        self.fig_ILS = plt.Figure(figsize=(6, 3), dpi=100)
        self.ax_ILS = self.fig_ILS.subplots(1, 1)
        self.ax_ILS.set_ylabel('Normalised intensity')
        self.ax_ILS.set_xlabel('Relative wavelength [nm]')
        self.ax_ILS.set_ylim([0, 1])
        self.ax_ILS.set_xlim([0, self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave])
        self.ax_ILS.grid(True)
        self.ax_ILS.plot([0, self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave], [0, 0], 'g', linewidth=1)
        self.ax_ILS.legend(('Instrument line shape',), loc=1, framealpha=1)
        self.fig_ILS.tight_layout()

        # Organise and draw canvas
        self.canv_ILS = FigureCanvasTkAgg(self.fig_ILS, master=self.frame_ILS)
        self.canv_ILS.draw()
        self.canv_ILS.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def dark_capture(self):
        """Controls dark spectrum capture"""
        if not self.__check_connection__():
            return

        # Set integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Get Spectrum
        self.wavelengths = self.spec_ctrl.wavelengths
        self.dark_spec = self.spec_ctrl.get_spec()

        # Correct calibation spectrum if we have one
        if self.cal_spec_raw is not None:
            self.cal_spec_corr = self.cal_spec_raw - self.dark_spec
            self.cal_spec_corr[self.cal_spec_corr < 0] = 0
            self.extract_ILS()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        self.dark_filename = '{}_dark.txt'.format(time)
        self.save_dark(self.save_path + self.dark_filename)

        # Change GUI label for dark file
        self.dark_file_label.configure(text=self.dark_filename)

        # Update plot with new data
        self.update_dark()

    def cal_capture(self):
        """Controls clear spectrum capture"""
        if not self.__check_connection__():
            return

        # Set integration time
        self.spec_ctrl.int_time = self.int_time.get()

        # Ignore first spectrum as it could have been acquired prior to
        self.cal_spec_raw = self.spec_ctrl.get_spec()

        # Dark subtract if we can
        if self.dark_spec is not None:
            self.cal_spec_corr = self.cal_spec_raw - self.dark_spec
            self.cal_spec_corr[self.cal_spec_corr < 0] = 0
            self.extract_ILS()

        # Generate filename based on time, then save file as text
        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        self.cal_filename = '{}_cal.txt'.format(time)
        self.save_cal_raw(self.save_path + self.cal_filename)

        # Change GUI label for dark file
        self.cal_file_label.configure(text=self.cal_filename)

        # Update plot with new data
        self.update_cal()

    def save_cal_raw(self, filename):
        """Save calibation data"""
        if self.wavelengths is None or self.cal_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.cal_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.cal_spec_raw]),
                   header='Raw calibration spectrum\n'
                          'Not dark-subtracted\n'
                          'Integration time={} ms\nWavelength [nm]\tIntensity [DN]'.format(self.int_time.get()))

    def save_dark(self, filename):
        """Save calibation data"""
        if self.wavelengths is None or self.dark_spec is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.wavelengths) != len(self.dark_spec):
            raise ValueError('Arrays are not the same length. Cannot save.')

        np.savetxt(filename, np.transpose([self.wavelengths, self.dark_spec]),
                   header='Dark spectrum\n'
                          'Integration time={} ms\nWavelength [nm]\tIntensity [DN]'.format(self.int_time.get()))

    def __check_connection__(self):
        """Checks spectrometer connection"""
        if isinstance(self.spec_ctrl, type(None)):
            print('No spectrometer found')
            return 0
        if self.spec_ctrl.spec is None:
            try:
                self.spec_ctrl.find_device()
            except SpectrometerConnectionError:
                print('No spectrometer found')
                return 0

        return 1    # Returns if connections were successful

    def update_dark(self):
        """Update dark plot with new data"""
        self.ax.lines[0].set_data(self.wavelengths, self.dark_spec)
        self.ax.set_xlim([self.wavelengths[0], self.wavelengths[-1]])
        self.canv.draw()

    def update_cal(self):
        """Update clear plot with new data"""
        self.ax.lines[1].set_data(self.wavelengths, self.cal_spec_raw)
        self.ax.set_xlim([self.wavelengths[0], self.wavelengths[-1]])
        self.canv.draw()

    def update_ILS_start(self):
        """Updates location of ILS extraction line, and performs extraction of ILS if calibration spectrum is present"""
        ILS_start = self.ILS_start.get()

        # Ensure the start of the stray range doesn't become less than the start
        if ILS_start >= self.ILS_end.get():
            ILS_start = self.ILS_end.get() - 0.1
            self.ILS_start.set(ILS_start)

        # Update plot
        self.min_ILS_line[0].set_data([ILS_start, ILS_start],
                                      [0, self.max_DN])
        self.ILS_range.set_x(ILS_start)
        self.ILS_range.set_width(self.ILS_end.get() - ILS_start)
        self.canv.draw()

        if self.cal_spec_corr is not None:
            self.extract_ILS()

    def update_ILS_end(self):
        """Updates location of ILS extraction line, and performs extraction of ILS if calibration spectrum is present"""
        ILS_end = self.ILS_end.get()

        # Ensure the end of the stray range doesn't become less than the start
        if ILS_end <= self.ILS_start.get():
            ILS_end = self.ILS_start.get() + 0.1
            self.ILS_end.set(ILS_end)

        # Update plot
        self.max_ILS_line[0].set_data([ILS_end, ILS_end],
                                        [0, self.max_DN])
        self.ILS_range.set_width(ILS_end - self.ILS_start.get())
        self.canv.draw()

        if self.cal_spec_corr is not None:
            self.extract_ILS()

    def choose_ILS(self):
        """Bring up filedialog to select ILS, then instigates loading"""
        # Bring up dialog to find file
        self.ILS_path = filedialog.askopenfilename(initialdir=self.save_path,
                                                   title='Select ILS file',
                                                   filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        self.load_ILS()

    def load_ILS(self):
        """Loads ILS from text file"""
        if not self.ILS_path:
            return

        # Update label in widget
        ILS_filename = self.ILS_path.split('/')[-1]
        if len(ILS_filename) > self.str_len_max:
            self.ILS_load_label.configure(text='...' + ILS_filename[-(self.str_len_max - 3):])
        else:
            self.ILS_load_label.configure(text=ILS_filename)

        # Extract data
        data = np.loadtxt(self.ILS_path)
        self.ILS_wavelengths, self.ILS = data.T

        # Update ILS plot
        self.plot_ILS()

    def save_ILS(self):
        """Saves ILS to text file"""
        if self.ILS is None or self.ILS_wavelengths is None:
            print('No instrument line shape data to save')
            return

        time = datetime.now().strftime('%Y-%m-%dT%H%M%S')
        ILS_filename = '{}_ILS.txt'.format(time)

        # Save file to numpy array
        np.savetxt(self.save_path + ILS_filename, np.transpose([self.ILS_wavelengths, self.ILS]),
                   header='Instrument line shape\n'
                          'Line start [nm]: {}\tLine End [nm]: {}\n'
                          'Wavelength [nm]\tIntensity [DN]'.format(self.ILS_start.get(), self.ILS_end.get()))

        # Update widget label
        self.ILS_save_label.configure(text=ILS_filename)

    def extract_ILS(self):
        """Extracts instrument line shape from the calibration spectrum, then draws it"""
        # Determine indices for extraction of ILS
        start_idx = np.argmin(np.absolute(self.wavelengths - self.ILS_start.get()))
        end_idx = np.argmin(np.absolute(self.wavelengths - self.ILS_end.get())) + 1

        # Extract ILS and normalise it
        # Need to use copy() here, otherwise we end up modifying cal_spec_corr when modifying ILS
        self.ILS = np.copy(self.cal_spec_corr[start_idx:end_idx])
        self.ILS /= np.amax(self.ILS)

        # Extract wavelengths then set to start at 0 nm
        self.ILS_wavelengths = self.wavelengths[start_idx:end_idx] - self.wavelengths[start_idx]

        # Update ILS plot
        self.plot_ILS()

    def plot_ILS(self):
        # Update plot
        self.ax_ILS.lines[0].set_data(self.ILS_wavelengths, self.ILS)
        self.ax_ILS.set_xlim([0, self.ILS_wavelengths[-1]])
        self.canv_ILS.draw()

        # Updates doas worker with ILS every time the plot is updated, so that we are in sync
        self.doas_worker.ILS = self.ILS

        # As we have a new ILS, we need to let the doas_worker that it has not been convolved with the ref spectrum
        self.doas_worker.ref_convolved = False


class RefPlot:
    """
    Plots and allows interaction with the reference spectrum
    """
    def __init__(self, frame, doas_worker=DOASWorker(), init_dir='C:\\', fig_size=(10,4), dpi=100, ref_spec_path=None):
        self.doas_worker = doas_worker

        self.init_dir = init_dir

        self.fig_size = fig_size
        self.dpi = dpi

        self.setts = SettingsGUI()
        self.pdx = 5
        self.pdy = 5

        # Setup gui
        self.__setup_gui__(frame)

        # If we have a pre-loaded reference spectrum plot it
        self.ref_spec_path = ref_spec_path
        if self.ref_spec_path is not None:
            self.load_ref_spec()

    def __setup_gui__(self, frame):
        # -------------------------
        # Reference Spectrum Setup
        # -------------------------
        self.frame = tk.LabelFrame(frame, text='Reference Spectrum', font=self.setts.mainFontBold,
                                      relief=tk.RAISED, borderwidth=2, bg=self.setts.bgColour)
        # self.refFrame.grid(row=0, column=0, padx=self.pdx, pady=self.pdy)
        self.frame.pack(side=tk.LEFT)

        self.loadRefFrame = tk.Frame(self.frame)
        self.loadRefFrame.grid(row=0, column=0, sticky='w')
        label = tk.Label(self.loadRefFrame, text='Filename:', font=self.setts.mainFontBold)
        label.grid(row=0, column=0, padx=self.pdx, pady=self.pdy)
        self.nameRef = tk.Label(self.loadRefFrame, text='None Selected', font=self.setts.mainFont)
        self.nameRef.grid(row=0, column=1, padx=self.pdx, pady=self.pdy)
        self.selectRef = ttk.Button(self.loadRefFrame, text='Load Spectrum', command=self.choose_ref_spec)
        self.selectRef.grid(row=0, column=2, padx=self.pdx, pady=self.pdy)

        self.conv_button = ttk.Button(self.frame, text='Convolve with ILS', command=self.conv_ref)
        self.conv_button.grid(row=0, column=1)

        # --------------------------------------------
        # FIGURE SETUP
        # --------------------------------------------
        self.FigRef = plt.Figure(figsize=self.fig_size, dpi=self.dpi)
        self.ax_SO2 = self.FigRef.add_subplot(111)

        self.ax_SO2.set_title('Reference spectrum: None')
        self.ax_SO2.set_ylabel(r'Absorption Cross Section [cm$^2$/molecule]')
        self.ax_SO2.set_xlabel('Wavelength [nm]')
        self.ax_SO2.tick_params(axis='both', direction='in', top='on', right='on')


        for i in range(2):
            self.ax_SO2.plot([250, 400], [0, 0], linewidth=1)
        self.ax_SO2.legend(('Reference', 'ILS-convolved'), loc=1, framealpha=1)

        self.refCanvas = FigureCanvasTkAgg(self.FigRef, master=self.frame)
        self.refCanvas.draw()
        self.refCanvas.get_tk_widget().grid(row=1, column=0, columnspan=2)
        # --------------------------------------------------------------------------------------------------------------

    def choose_ref_spec(self):
        """Load reference spectrum"""
        self.ref_spec_path = filedialog.askopenfilename(initialdir=self.init_dir,
                                                        title='Select reference spectrum',
                                                        filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        self.load_ref_spec()

    def load_ref_spec(self):
        """Loads reference spectrum"""
        if not self.ref_spec_path:
            return
        if len(self.ref_spec_path) > 53:
            self.nameRef.configure(text='...' + self.ref_spec_path[-50:])
        else:
            self.nameRef.configure(text=self.ref_spec_path)
        self.doas_worker.load_ref_spec(self.ref_spec_path, 'SO2')

        self.ref_spec_file = self.ref_spec_path.split('/')[-1]

        # Set convolution plot to zero then update
        self.ax_SO2.lines[1].set_data([self.doas_worker.ref_spec['SO2'][0, 0], self.doas_worker.ref_spec['SO2'][-1,0]],
                                      [0,0])
        # Plot reference spectrum
        self.plot_ref_spec()

    def plot_ref_spec(self):
        """Plot up reference spectrum"""
        self.ax_SO2.lines[0].set_data(self.doas_worker.ref_spec['SO2'][:, 0], self.doas_worker.ref_spec['SO2'][:, 1])
        self.ax_SO2.set_xlim([self.doas_worker.ref_spec['SO2'][0, 0], self.doas_worker.ref_spec['SO2'][-1, 0]])
        self.ax_SO2.set_ylim([0, np.amax(self.doas_worker.ref_spec['SO2'][:, 1])])
        self.ax_SO2.set_title('Reference Spectrum: %s' % self.ref_spec_file)
        self.refCanvas.draw()

    def conv_ref(self):
        if self.doas_worker.ILS is None:
            print('No ILS to convolved reference spectrum with')
            return

        # If reference spectrum dictionary is empty - return
        if not self.doas_worker.ref_spec:
            print('No reference spectrum loaded')
            return

        # Convolve reference spectrum
        self.doas_worker.conv_ref_spec()

        # Update plot (only if convolution was successful - it fails if wavelength data is not present - need to acquire spectrum with spectrometer first)
        if 'SO2' in self.doas_worker.ref_spec_conv:
            self.ax_SO2.lines[1].set_data(self.doas_worker.wavelengths, self.doas_worker.ref_spec_conv['SO2'])
            self.refCanvas.draw()