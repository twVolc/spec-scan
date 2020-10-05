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
from save_spec import SaveSpectra

import numpy as np
from datetime import datetime
import os
import time
from PIL import Image, ImageTk

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
        self.scan_cont = ScanProperties()
        self.save_obj = SaveSpectra(self.doas_worker, self.scan_proc, None)

        self.str_len_max = 30
        self.pdx = 5
        self.pdy = 5

        self.init_dir = init_dir    # Directory to begin filedialog prompts
        self.save_path = init_dir   # Directory to save processed date to
        self.dark_path = None
        self.clear_path = None
        self.plume_path = None
        self.scan_files = None
        self.scan_dir = None
        self.doas_path = None

        self.scan_start = True   # Used to adjust button to control whether new scan is being processed or whether stepping to next spec

        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):
        self.frame = tk.LabelFrame(frame, text='Post-Processing', relief=tk.RAISED, borderwidth=5)

        # Dark load
        label = ttk.Label(self.frame, text='Dark spectrum file:')
        self.dark_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_dark_butt = ttk.Button(self.frame, text='Load Dark', command=self.select_dark)

        row = 0
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.dark_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_dark_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')

        # Clear load
        label = ttk.Label(self.frame, text='Clear spectrum file:')
        self.clear_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_clear_butt = ttk.Button(self.frame, text='Load Clear', command=self.select_clear)

        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.clear_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_clear_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')

        # Plume load
        label = ttk.Label(self.frame, text='Plume spectrum file:')
        self.plume_label = ttk.Label(self.frame, text='N/A', width=25)
        self.load_plume_butt = ttk.Button(self.frame, text='Load Plume', command=self.select_plume)

        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.plume_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_plume_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')
        row += 1


        self.scan_frame = ttk.LabelFrame(self.frame, text='Scan', relief=tk.RAISED)
        self.scan_frame.grid(row=row, column=0, columnspan=2, sticky='nsew')

        row += 1
        self.process_butt = ttk.Button(self.frame, text='Process', command=self.process_doas)
        self.process_butt.grid(row=row, column=0, columnspan=2, sticky='nsew', padx=self.pdx, pady=self.pdy)

        # Load scan
        row=0


        self.batch_proc = tk.IntVar()
        self.batch_proc.set(1)
        self.batch_check = ttk.Checkbutton(self.scan_frame, text='Batch process', variable=self.batch_proc)
        self.batch_check.grid(row=row, column=0, padx=self.pdx, pady=self.pdy, sticky='w')
        row += 1

        label = ttk.Label(self.scan_frame, text='Scan directory:')
        self.scan_label = ttk.Label(self.scan_frame, text='N/A', width=25)
        self.load_scan_butt = ttk.Button(self.scan_frame, text='Load Scan', command=self.load_scan)

        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.scan_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)
        row += 1
        self.load_scan_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')
        row += 1

        label = ttk.Label(self.scan_frame, text='Select clear spectrum:')
        self.clear_spec_var = tk.StringVar()
        self.clear_spec_select = ttk.Combobox(self.scan_frame, width=27, justify='center',
                                              textvariable=self.clear_spec_var, state='readonly')
        self.clear_spec_select['values'] = ["--load scan directory--"]
        self.clear_spec_select.current(0)
        self.clear_spec_select.bind("<<ComboboxSelected>>", self.update_clear_combo)

        label.grid(row=row, column=0, padx=self.pdx, pady=self.pdy)
        self.clear_spec_select.grid(row=row, column=1, padx=self.pdx, pady=self.pdy)



    def select_dark(self):
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

        self.load_dark()

    def load_dark(self):
        """Load dark spectrum and update plots"""
        # Extract data
        data = np.loadtxt(self.dark_path)
        self.doas_worker.wavelengths, self.doas_worker.dark_spec = data.T

        # Update dark plot
        self.spec_plot.update_dark()

    def select_clear(self):
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

        # Update the clear spectrum doas_worker + plot
        self.load_clear()

    def load_clear(self):
        """Updates clear spectrum by loading data and updating plots (used by a few functions)"""
        # Extract data
        if os.path.exists(self.clear_path):
            data = np.loadtxt(self.clear_path)
            self.doas_worker.wavelengths, self.doas_worker.clear_spec_raw = data.T

            # Update dark plot
            self.spec_plot.update_clear()
        else:
            print('The file could not be loaded, please select a file that exists')

    def update_clear_combo(self, event):
        """updates clear spectrum from combobox"""
        # Update clear spectrum
        try:
            self.clear_path = self.scan_dir + self.clear_spec_var.get()
        except TypeError:
            return

        # Check path exists, as this function is called even if using combobox before a scan directory has been loaded
        if os.path.exists(self.clear_path):
            self.load_clear()
        else:
            return

    def select_plume(self):
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

        self.load_plume()

    def load_plume(self, save=True):
        """Loads in plume spectrum and updates plot. Also tries to process the data"""
        # Extract data
        data = np.loadtxt(self.plume_path)
        self.doas_worker.wavelengths, self.doas_worker.plume_spec_raw = data.T

        # Update dark plot
        self.spec_plot.update_plume()

        # Try processing data
        self.doas_worker.process_doas()

        # Update plot and save processed spectrum if data was processed
        if self.doas_worker.processed_data:
            self.doas_plot.update_plot()
            if save:
                self.set_doas_filename(self.plume_path.rsplit('/')[-1])
                self.save_obj.save_processed_spec(self.dark_path, self.clear_path, self.plume_path, self.doas_path)


    def load_scan(self):
        """Load scan directory"""
        # Save old scan directory if we need to revert back
        scan_dir_old = self.scan_dir

        # Open dialog to select new scan directory
        self.scan_dir = filedialog.askdirectory(initialdir=self.init_dir, title='Select scan folder')

        # Revert to old scan directory if file dialog box was closed without selecting a directory
        if not self.scan_dir:
            self.scan_dir = scan_dir_old
            return

        # Configure scan directory label widget - truncate string if it is too long for the widget
        self.scan_dir += '/'
        self.init_dir = self.scan_dir
        if len(self.scan_dir) > self.str_len_max:
            self.scan_label.configure(text='...' + self.scan_dir[-(self.str_len_max - 8):])
        else:
            self.scan_label.configure(text=self.scan_dir)

        all_files = os.listdir(self.scan_dir)

        dark_files = [f for f in all_files if 'dark.txt' in f.lower()]
        self.scan_files = [f for f in all_files if 'plume.txt' in f.lower() or 'clear.txt' in f.lower()]     # Includes any clear files in the directory

        # For now just work with one dark spectrum
        if dark_files:
            self.dark_path = self.scan_dir + dark_files[0]
            self.load_dark()

        if self.scan_files:
            # Organise clear spectrum combobox
            self.clear_spec_select['values'] = self.scan_files
            self.clear_spec_select.current(0)
            self.update_clear_combo(0)

    def __update_scan_params__(self):
        """Updates scan paramters"""
        self.scan_proc.plume_speed = self.cd_plot.plume_speed
        self.scan_proc.plume_distance = self.cd_plot.plume_dist

    def process_doas(self):
        """Perform doas retrieval on loaded data"""
        # No need to check plume_path, as clear_path comes from the first plume spectrum in scan_dir (this is not the case if loading individual spectra)
        if self.dark_path is None or not self.scan_files:
            print('Scan directory must contain dark and plume spectra to be processed')
            return

        # Create directory to save processed data
        self.save_path = self.scan_dir + 'Processing_1/'
        i = 2
        while os.path.exists(self.save_path):
            self.save_path = self.scan_dir + 'Processing_{}/'.format(i)
            i += 1
        os.mkdir(self.save_path)

        # Initialise scan_proc paramters
        self.__update_scan_params__()   # Gets plume speed and distance
        self.scan_proc.scan_angles = np.array([])
        self.scan_proc.column_densities = np.array([])
        scan_angle = 0

        self.cd_plot.clear_plot()

        # Update clear again, just in case other clear files have been loaded
        self.update_clear_combo(0)

        batch = self.batch_proc.get()
        if batch:
            save = True
        else:
            save = False

        # Process all automatically
        for plume_file in self.scan_files:
            # Load plume spectrum and automatically process it (processing is handled outside of this class)
            self.plume_path = self.scan_dir + plume_file
            self.load_plume(save=save)

            if not batch:
                # Wait for button click to proceed to next spectrum
                var = tk.IntVar()
                self.process_butt.configure(text='NEXT SPECTRUM', command=lambda: var.set(1))
                self.process_butt.wait_variable(var)

            # Update CD scan plot
            scan_angle += self.scan_cont.scan_incr
            self.scan_proc.add_data(scan_angle, self.doas_worker.column_density['SO2'])
            self.cd_plot.update_plot()

            if not batch:
                # Save processed spectrum
                self.set_doas_filename(plume_file)
                self.save_obj.save_processed_spec(self.dark_path, self.clear_path, self.plume_path, self.doas_path)

        # Calculate emission rate
        self.scan_proc.calc_emission_rate()
        self.cd_plot.update_emission_rate()

        # Save scan data
        self.save_obj.save_scan(self.save_path)

        if not batch:
            # Re configure the button to control the start of a processing loop again
            self.process_butt.configure(text='Process', command=self.process_doas)

            # Add a small lag when processing has finished
            # Prevents the user quickly moving to the next scan by accident
            time.sleep(1)

    def set_doas_filename(self, filename):
        """Organise doas filename for saving purposes"""
        filename = filename.split('.')[0].replace('clear', 'doas')
        filename = filename.split('.')[0].replace('plume', 'doas')  # Filename may contain plume or clear
        doas_filename = filename + '_1.txt'
        self.doas_path = self.save_path + doas_filename

        idx = 2
        while os.path.exists(self.doas_path):
            self.doas_path = self.save_path + filename + '_{}.txt'.format(idx)
            idx += 1

        # try:
        #     np.savetxt(self.doas_path, np.transpose([self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit['SO2'],
        #                                              self.doas_worker.abs_spec_cut]),
        #                header='Processed DOAS spectrum\n'
        #                       'Dark spectrum: {}\nClear spectrum: {}\nPlume spectrum: {}\n'
        #                       'Shift: {}\nStretch: {}\n'
        #                       'Stray range [nm]: {}:{}\nFit window [nm]: {}:{}\n'
        #                       'Column density [ppm.m]: {}\n'
        #                       'Wavelength [nm]\tReference spectrum (fitted)\tAbsorbance spectrum'.format(
        #                       self.dark_path, self.clear_path, self.plume_path,
        #                       self.doas_worker.shift, self.doas_worker.stretch,
        #                       self.doas_worker.start_stray_wave, self.doas_worker.end_stray_wave,
        #                       self.doas_worker.start_fit_wave, self.doas_worker.end_fit_wave,
        #                       self.doas_worker.column_density))
        # except:
        #     print('Permission denied! Could not save processed spectra. Please change the save directory')

    def save_scan(self):
        """Saves scan information"""
        filename = 'Scan_data.txt'

        try:
            np.savetxt(self.save_path + filename, np.transpose([self.scan_proc.scan_angles,
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


class DirectoryWatcherFrame:
    """
    Class to create widget for watchign a directory and analysing data as it arrives
    """
    def __init__(self, frame, doas_worker=DOASWorker(), init_dir='C:\\', generate_frame=True):
        self.frame = frame
        self.doas_worker = doas_worker
        self.init_dir = init_dir
        self.watching = False

        self.pdx = 5
        self.pdy = 5
        self.str_len_max = 30

        self.img_size = (40, 40)
        self.img_on = ImageTk.PhotoImage(Image.open('./icons/green-led.png').resize(self.img_size, Image.ANTIALIAS))
        self.img_off = ImageTk.PhotoImage(Image.open('./icons/red-led.png').resize(self.img_size, Image.ANTIALIAS))

        if generate_frame:
            self.generate_frame()

    def generate_frame(self):
        """Build gui frame"""
        self.frame = tk.LabelFrame(self.frame, text='Directory watcher', relief=tk.RAISED, borderwidth=5)

        row = 0

        # Define whether to get plume parameters from file
        self._auto_plume_params = tk.IntVar()
        check = ttk.Checkbutton(self.frame, text='Auto-read plume parameters', variable=self._auto_plume_params,
                                command=self.update_param)
        check.grid(row=row, column=0, columnspan=2, sticky='w', padx=self.pdx, pady=self.pdy)

        label = ttk.Label(self.frame, text='Watch directory:')
        self.watch_label = ttk.Label(self.frame, text='N/A', width=25)
        self.select_butt = ttk.Button(self.frame, text='Select directory', command=self.select_dir)
        row += 1
        label.grid(row=row, column=0, sticky='e', padx=self.pdx, pady=self.pdy)
        self.watch_label.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')
        row += 1
        self.select_butt.grid(row=row, column=1, padx=self.pdx, pady=self.pdy, sticky='nsew')
        self.indicator = tk.Canvas(self.frame, width=self.img_size[0], height=self.img_size[1],
                                   bd=0, highlightthickness=0)
        self.indicator.create_image(0, 0, image=self.img_off, anchor='nw', tags='IMG')
        self.indicator.grid(row=row, column=0, rowspan=2, sticky='e')
        row += 2

        self.frame.grid_columnconfigure(1, weight=1)

        self.watch_butt = ttk.Button(self.frame, text='Start Watching', command=self.watch_directory)
        self.watch_butt.grid(row=row, column=0, columnspan=2, sticky='nsew', padx=self.pdx, pady=self.pdy)

    def select_dir(self):
        """Selects directory to be watched"""
        watch_path = filedialog.askdirectory(initialdir=self.init_dir, title='Select watch directory')
        if not watch_path:
            return

        self.init_dir = watch_path

        # Configure watch path label
        if len(watch_path) > self.str_len_max:
            self.watch_label.configure(text='...' + watch_path[-(self.str_len_max - 3):])
        else:
            self.watch_label.configure(text=watch_path)

        # Update doas worker to new path
        self.doas_worker.watch_dir = watch_path

    def watch_directory(self):
        """Initiates directory watcher to begin processing anythign which enters the directory"""
        if not self.watching:
            self.doas_worker.start_continuous_processing()
            self.watch_butt.configure(text='Stop Watching')
            self.watching = True
            self.indicator.delete('IMG')
            self.indicator.create_image(0, 0, image=self.img_on, anchor='nw', tags='IMG')


        elif self.watching:
            self.doas_worker.stop_continuous_processing()
            self.watch_butt.configure(text='Start Watching')
            self.watching = False
            self.indicator.delete('IMG')
            self.indicator.create_image(0, 0, image=self.img_off, anchor='nw', tags='IMG')

        else:
            raise AttributeError('Encountered an unexpected value for watching bool')

    def update_param(self):
        """Updates doas_worker parameter for auto plume parameter loading from a text file"""
        self.doas_worker.auto_plume_params = self.auto_plume_params

    @property
    def auto_plume_params(self):
        return bool(self._auto_plume_params.get())
