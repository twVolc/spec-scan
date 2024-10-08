import tkinter as tk
import tkinter.ttk as ttk

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.dates import date2num, DateFormatter, HourLocator

import numpy as np
import queue
import datetime

from gui_subs import SettingsGUI
from doas_routine import DOASWorker, ScanProcess
from acquisition_gui import AcquisitionFrame

plt.style.use('dark_background')

refresh_rate = 100

class SpectraPlot:
    """
    Generates a widget containing 3 subplots of spectra -> dark, clear (Fraunhofer), in-plume
    """

    def __init__(self, root, frame, doas_worker=DOASWorker(), doas_plot=None, figsize=(10, 3), dpi=100):
        self.root = root
        self.setts = SettingsGUI()
        self.doas_worker = doas_worker
        self.doas_plot = doas_plot
        self.doas_worker.fig_spec = self

        self.figsize = figsize
        self.dpi = dpi

        self.max_DN = 2**16 - 1  # Maximum DN for spectrometer

        # Could use threads and queues to update plots, or just us simple functions which Acquisition frame calls
        self.Q = queue.Queue()

        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):
        """Organise widget"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=4)

        # -----------------------------------------------
        # WIDGET SETUP
        # ------------------------------------------------
        self.frame2 = ttk.Frame(self.frame)
        self.frame2.pack(side=tk.TOP, fill=tk.X, expand=1)

        # STRAY RANGE
        self.stray_start = tk.DoubleVar()
        self.stray_start.set(self.doas_worker.start_stray_wave)
        self.stray_box_start = tk.Spinbox(self.frame2, from_=0, to=400, increment=0.1, width=5,
                                             textvariable=self.stray_start, command=self.update_stray_start)
        self.stray_end = tk.DoubleVar()
        self.stray_end.set(self.doas_worker.end_stray_wave)
        self.stray_box_end = tk.Spinbox(self.frame2, from_=1, to=400, increment=0.1, width=5,
                                           textvariable=self.stray_end, command=self.update_stray_end)

        label = tk.Label(self.frame2, text='Stray light correction (min.):').pack(side=tk.LEFT)
        self.stray_box_start.pack(side=tk.LEFT)
        label = tk.Label(self.frame2, text='Stray light correction (max.):').pack(side=tk.LEFT)
        self.stray_box_end.pack(side=tk.LEFT)

        # FIT WINDOW
        self.fit_wind_start = tk.DoubleVar()
        self.fit_wind_start.set(self.doas_worker.start_fit_wave)
        self.fit_wind_box_start = tk.Spinbox(self.frame2, from_=0, to=400, increment=0.1, width=5,
                                             textvariable=self.fit_wind_start, command=self.update_fit_wind_start)
        self.fit_wind_end = tk.DoubleVar()
        self.fit_wind_end.set(self.doas_worker.end_fit_wave)
        self.fit_wind_box_end = tk.Spinbox(self.frame2, from_=1, to=400, increment=0.1, width=5,
                                           textvariable=self.fit_wind_end, command=self.update_fit_wind_end)

        self.fit_wind_box_end.pack(side=tk.RIGHT)
        label = tk.Label(self.frame2, text='Fit wavelength (max.):').pack(side=tk.RIGHT)
        self.fit_wind_box_start.pack(side=tk.RIGHT)
        label = tk.Label(self.frame2, text='Fit wavelength (min.):').pack(side=tk.RIGHT)

        # ------------------------------------------------
        # FIGURE SETUP
        # ------------------------------------------------
        self.fig = plt.Figure(figsize=self.figsize, dpi=self.dpi)

        self.ax = self.fig.subplots(1, 1)
        self.ax.set_ylabel('DN')
        self.ax.set_ylim([0, self.max_DN])
        self.ax.set_xlim([250, 400])
        self.ax.set_xlabel('Wavelength [nm]')
        self.ax.grid(True)
        self.plt_colours = ['b', 'g', 'r']
        for i in range(3):
            self.ax.plot([250, 400], [0, 0], self.plt_colours[i], linewidth=1)
        self.ax.legend(('Dark', 'Clear', 'Plume'), loc=2, framealpha=1)
        self.fig.tight_layout()

        # Stray light
        self.min_stray_line = self.ax.plot([self.doas_worker.start_stray_wave,self.doas_worker.start_stray_wave],
                                     [0, self.max_DN], 'm')
        self.max_stray_line = self.ax.plot([self.doas_worker.end_stray_wave,self.doas_worker.end_stray_wave],
                                     [0, self.max_DN], 'm')

        width = self.doas_worker.end_stray_wave - self.doas_worker.start_stray_wave
        self.stray_range = Rectangle((self.doas_worker.start_stray_wave,0), width=width, height=self.max_DN,
                                  fill=True, facecolor='m', alpha=0.2)
        self.stray_range_patch = self.ax.add_patch(self.stray_range)

        # Fit window
        self.min_line = self.ax.plot([self.doas_worker.start_fit_wave, self.doas_worker.start_fit_wave],
                                     [0, self.max_DN], 'y')
        self.max_line = self.ax.plot([self.doas_worker.end_fit_wave, self.doas_worker.end_fit_wave],
                                     [0, self.max_DN], 'y')

        width = self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave
        self.fit_wind = Rectangle((self.doas_worker.start_fit_wave, 0), width=width, height=self.max_DN,
                                  fill=True, facecolor='y', alpha=0.2)
        self.fit_wind_patch = self.ax.add_patch(self.fit_wind)

        # Organise and draw canvas
        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Instigate canvas drawing worker
        self.__draw_canv__()

    def update_dark(self):
        """Update dark plot with new data"""
        self.ax.lines[0].set_data(self.doas_worker.wavelengths, self.doas_worker.dark_spec)
        self.ax.set_xlim([self.doas_worker.wavelengths[0], self.doas_worker.wavelengths[-1]])
        self.Q.put(1)

    def update_clear(self):
        """Update clear plot with new data"""
        self.ax.lines[1].set_data(self.doas_worker.wavelengths, self.doas_worker.clear_spec_raw)
        self.ax.set_xlim([self.doas_worker.wavelengths[0], self.doas_worker.wavelengths[-1]])
        self.Q.put(1)

    def update_plume(self):
        """Update clear plot with new data"""
        self.ax.lines[2].set_data(self.doas_worker.wavelengths, self.doas_worker.plume_spec_raw)
        self.ax.set_xlim([self.doas_worker.wavelengths[0], self.doas_worker.wavelengths[-1]])
        self.Q.put(1)

    def update_stray_start(self):
        """Updates stray light range on plot"""
        stray_start = self.stray_start.get()

        # Ensure the end of the stray range doesn't become less than the start
        if stray_start >= self.stray_end.get():
            stray_start = self.stray_end.get() - 0.1
            self.stray_start.set(stray_start)

        # Update DOASWorker with new stray range
        self.doas_worker.start_stray_wave = stray_start

        # Update plot
        self.min_stray_line[0].set_data([self.doas_worker.start_stray_wave, self.doas_worker.start_stray_wave], [0, self.max_DN])
        self.stray_range.set_x(self.doas_worker.start_stray_wave)
        self.stray_range.set_width(self.doas_worker.end_stray_wave - self.doas_worker.start_stray_wave)
        self.Q.put(1)

        # Process doas if we can (if it has previously been processed so we know we have all the data we need)
        if self.doas_worker.processed_data:
            self.doas_worker.stray_corrected = False
            self.doas_worker.process_doas()
            self.doas_plot.update_plot()

    def update_stray_end(self):
        """Updates stray light range on plot"""
        stray_end = self.stray_end.get()

        # Ensure the end of the stray range doesn't become less than the start
        if stray_end <= self.stray_start.get():
            stray_end = self.stray_start.get() + 0.1
            self.stray_end.set(stray_end)

        # Update DOASWorker with new stray range
        self.doas_worker.end_stray_wave = stray_end

        # Update plot
        self.max_stray_line[0].set_data([self.doas_worker.end_stray_wave, self.doas_worker.end_stray_wave], [0, self.max_DN])
        self.stray_range.set_width(self.doas_worker.end_stray_wave - self.doas_worker.start_stray_wave)
        self.Q.put(1)

        if self.doas_worker.processed_data:
            self.doas_worker.stray_corrected = False
            self.doas_worker.process_doas()
            self.doas_plot.update_plot()

    def update_fit_wind_start(self):
        """updates fit window on plot"""
        fit_wind_start = self.fit_wind_start.get()

        # Ensure the end of the fit window doesn't become less than the start
        if fit_wind_start >= self.fit_wind_end.get():
            fit_wind_start = self.fit_wind_end.get() - 0.1
            self.fit_wind_start.set(fit_wind_start)

        # Update DOASWorker with new fit window
        self.doas_worker.start_fit_wave = fit_wind_start

        # Update plot
        self.min_line[0].set_data([self.doas_worker.start_fit_wave, self.doas_worker.start_fit_wave], [0, self.max_DN])
        self.fit_wind.set_x(self.doas_worker.start_fit_wave)
        self.fit_wind.set_width(self.doas_worker.end_fit_wave-self.doas_worker.start_fit_wave)
        self.Q.put(1)

        if self.doas_worker.processed_data:
            self.doas_worker.process_doas()
            self.doas_plot.update_plot()

    def update_fit_wind_end(self):
        """updates fit window on plot"""
        fit_wind_end = self.fit_wind_end.get()

        # Ensure the end of the fit window doesn't become less than the start
        if fit_wind_end <= self.fit_wind_start.get():
            fit_wind_end = self.fit_wind_start.get() + 0.1
            self.fit_wind_end.set(fit_wind_end)

        # Update DOASWorker with new fit window
        self.doas_worker.end_fit_wave = fit_wind_end

        # Update plot
        self.max_line[0].set_data([self.doas_worker.end_fit_wave, self.doas_worker.end_fit_wave], [0, self.max_DN])
        self.fit_wind.set_width(self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave)
        self.Q.put(1)

        if self.doas_worker.processed_data:
            self.doas_worker.process_doas()
            self.doas_plot.update_plot()

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            # print('Got {} from Q'.format(update))
            if update == 1:
                self.canv.draw()
            else:
                print('Closing canvas drawing')
                return
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

    def close_widget(self):
        """Closes widget cleanly, by stopping __draw_canv__()"""
        self.Q.put(2)
        # print('Added to Q')


class DOASPlot:
    """
    Generates a widget containing the DOAS fit plot
    """
    def __init__(self, root, frame, doas_worker=DOASWorker(), figsize=(10, 3), dpi=100, species='SO2'):
        self.root = root

        self.setts = SettingsGUI()
        self.doas_worker = doas_worker
        self.doas_worker.fig_doas = self

        self.species = species

        self.figsize = figsize
        self.dpi = dpi

        self.Q = queue.Queue()

        self.acq_obj = None

        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):
        """Organise widget"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=4)

        # -----------------------------------------------
        # WIDGET SETUP
        # ------------------------------------------------
        self.frame2 = ttk.Frame(self.frame)
        self.frame2.pack(side=tk.TOP, fill=tk.X, expand=1)

        label = tk.Label(self.frame2, text='Shift spectrum:').pack(side=tk.LEFT)
        # label.grid(row=0, column=0)
        self.shift = tk.IntVar()
        self.shift.set(self.doas_worker.shift)
        self.shift_box = tk.Spinbox(self.frame2, from_=-20, to=20, increment=1, width=3,
                                             textvariable=self.shift, command=self.update_shift)
        # self.fit_wind_box_start.grid(row=0, column=1)
        self.shift_box.pack(side=tk.LEFT)

        label2 = tk.Label(self.frame2, text='Stretch spectrum:').pack(side=tk.LEFT)
        # label2.grid(row=0, column=2)
        self.stretch = tk.IntVar()
        self.stretch.set(self.doas_worker.stretch)
        self.stretch_box = tk.Spinbox(self.frame2, from_=-999, to=999, increment=1, width=4,
                                           textvariable=self.stretch, command=self.update_stretch)
        # self.fit_wind_box_end.grid(row=0, column=3)
        self.stretch_box.pack(side=tk.LEFT)

        # Save button
        self.save_butt = ttk.Button(self.frame2, text='Save spectra', command=self.save_spectra)
        self.save_butt.pack(side=tk.RIGHT, anchor='e')

        # ----------------------------------------------------------------
        # Tabbed figure setup for all species and residual
        # ----------------------------------------------------------------
        # Setup up tab wideget for each species
        self.tabs = ttk.Notebook(self.frame)
        self.tabs.bind('<Button-1>', self.__update_tab__)
        self.species_tabs = dict()
        self.species_tabs['Total'] = ttk.Frame(self.tabs, borderwidth=2)
        self.tabs.add(self.species_tabs['Total'], text='Total')
        for spec in self.species:
            self.species_tabs[spec] = ttk.Frame(self.tabs, borderwidth=2)
            self.tabs.add(self.species_tabs[spec], text=spec)
        self.species_tabs['residual'] = ttk.Frame(self.tabs, borderwidth=2)
        self.tabs.add(self.species_tabs['residual'], text='residual')
        self.tabs.pack(side=tk.TOP, fill="both", expand=1)

        # Generate DOASFigure object for each species
        self.species_plots = dict()
        for spec in self.species:
            self.species_plots[spec] = DOASFigure(self.species_tabs[spec], self.doas_worker, spec,
                                                  self.figsize, self.dpi)

        # Generate DOASFigure object for full absorbance
        self.species_plots['Total'] = DOASFigure(self.species_tabs['Total'], self.doas_worker, 'Total',
                                                    self.figsize, self.dpi)

        # Generate DOASFigure object for residual
        self.species_plots['residual'] = DOASFigure(self.species_tabs['residual'], self.doas_worker, 'residual',
                                                    self.figsize, self.dpi)
        # Instigate canvas drawing worker
        self.__draw_canv__()

    def update_shift(self):
        """Updates DOASWorker shift value for aligning spectra"""
        # Set shift in DOASWorker object
        self.doas_worker.shift = self.shift.get()

        # If we have a processed spectrum we now must update it
        if self.doas_worker.processed_data:
            self.doas_worker.process_doas()
            self.update_plot()

    def update_stretch(self):
        """Updates DOASWorker stretch value for aligning spectra"""
        self.doas_worker.stretch = self.stretch.get()

        # If we have a processed spectrum we now must update it
        if self.doas_worker.processed_data:
            self.doas_worker.process_doas()
            self.update_plot()

    def update_plot(self):
        """Updates doas plot"""
        for plot in self.species_plots:
            self.species_plots[plot].update_plot()

        # Draw updates
        self.Q.put(1)

    def __update_tab__(self, event):
        """
        Controls drawing of tab canvas when tab is selected
        Drawing in this manner reduces lag when playing with a figure, as only the current figure is drawn
        :return:
        """
        self.Q.put(1)

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            # print('Got message')
            if update == 1:
                # Only draw canvas of currently selected tab (save processing power)
                species = self.tabs.tab(self.tabs.select(), "text")
                self.species_plots[species].canv.draw()
            else:
                return
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

    def close_widget(self):
        """Closes widget cleanly, by stopping __draw_canv__()"""
        self.Q.put(2)

    def save_spectra(self):
        """Saves processed DOAS spectra"""
        if isinstance(self.acq_obj, AcquisitionFrame):
            print('Saving...')
            self.acq_obj.save_processed_spec()


class DOASFigure:
    """
    Class for generating a DOAS-style figure with absorbance and reference spectrum fitting
    """
    def __init__(self, frame, doas_worker, species, figsize, dpi):
        self.frame = frame
        self.doas_worker = doas_worker
        self.species = species
        self.figsize = figsize
        self.dpi = dpi
        # ------------------------------------------------
        # FIGURE SETUP
        # ------------------------------------------------
        self.fig = plt.Figure(figsize=self.figsize, dpi=self.dpi)

        self.ax = self.fig.subplots(1, 1)
        self.ax.set_ylabel('Absorbance')
        self.ax.set_ylim([-0.2, 0.2])
        self.ax.set_xlim([self.doas_worker.start_fit_wave, self.doas_worker.end_fit_wave])
        self.ax.set_xlabel('Wavelength [nm]')
        self.ax.grid(True)
        self.plt_colours = ['b', 'r']
        for i in range(2):
            self.ax.plot([250, 400], [0, 0], self.plt_colours[i], linewidth=1)
        self.ax.legend(('Measured', 'Fitted reference'), loc=1, framealpha=1)
        self.ax.set_title('SO2 Column density [ppm.m]: N/A          STD Error: N/A')
        self.fig.tight_layout()

        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def update_plot(self):
        """Updates doas plot"""
        # Update plot lines with new data
        # self.ax.lines[0].set_data(self.doas_worker.wavelengths_cut, self.doas_worker.abs_spec_cut[self.species])
        self.ax.lines[0].set_data(self.doas_worker.wavelengths_cut, self.doas_worker.abs_spec_species[self.species])

        if self.species != 'residual':
            self.ax.lines[1].set_data(self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit[self.species])

        # Set axis limits
        self.ax.set_xlim([self.doas_worker.wavelengths_cut[0], self.doas_worker.wavelengths_cut[-1]])
        ylims = np.amax(np.absolute(self.doas_worker.abs_spec_species[self.species]))
        ylims *= 1.15
        if ylims == 0:
            ylims = 0.05
        self.ax.set_ylim([-ylims, ylims])
        self.ax.set_title('SO2 Column density [ppm.m]: {}          STD Error: {}'.format(
            self.doas_worker.column_density['SO2'], self.doas_worker.std_err))


class CDPlot:
    """
    Class to plot column densities retrieved from a scan or traverse sequence
    """
    def __init__(self, root, frame, doas_worker=DOASWorker(), fig_size=(5,5), dpi=100):
        self.root = root

        self.doas_worker = doas_worker
        self.doas_worker.fig_scan = self
        self.scan_proc = self.doas_worker.scan_proc
        self.stds = []  # List of std lines

        self.Q = queue.Queue()

        self.fig_size=fig_size
        self.dpi = dpi

        self.x_ax_min = 20

        self.__setup_gui__(frame)

    def __setup_gui__(self,frame):
        """Controls widget setup"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=5)

        # INPUT FRAME SETUP
        self.frame_inputs = ttk.Frame(self.frame)
        self.frame_inputs.pack(side=tk.LEFT, anchor='nw')

        label_dist = ttk.Label(self.frame_inputs, text='Plume distance [m]:')
        self._plume_dist = tk.IntVar()
        self.plume_dist = 1000
        self.plume_dist_box = ttk.Entry(self.frame_inputs, width=5, textvariable=self._plume_dist)

        label_dist.grid(row=0, column=0, padx=5, pady=5, sticky='e')
        self.plume_dist_box.grid(row=0, column=1, padx=5, pady=5)

        label_speed = ttk.Label(self.frame_inputs, text='Plume speed [m/s]:')
        self._plume_speed = tk.DoubleVar()
        self.plume_speed = 5.0
        self.plume_speed_box = ttk.Entry(self.frame_inputs, width=5, textvariable=self._plume_speed)

        label_speed.grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.plume_speed_box.grid(row=1, column=1, padx=5, pady=5)

        label = ttk.Label(self.frame_inputs, text='Scan emission rate [kg/s]:').grid(
            row=2, column=0, padx=5, pady=5, sticky='e')
        self.emission_rate = ttk.Label(self.frame_inputs, text='N/A')
        self.emission_rate.grid(row=2, column=1, padx=5, pady=5)

        label = ttk.Label(self.frame_inputs, text='Scan emission rate [t/day]:').grid(
            row=3, column=0, padx=5, pady=5, sticky='e')
        self.emission_rate_td = ttk.Label(self.frame_inputs, text='N/A')
        self.emission_rate_td.grid(row=3, column=1, padx=5, pady=5)

        # FIGURE SETUP
        self.fig = plt.Figure(figsize=self.fig_size, dpi=self.dpi)
        self.ax = self.fig.subplots(1,1)
        self.ax.set_ylabel('Column Density [ppm.m]')
        self.ax.set_xlabel('Scan Angle [deg]')
        self.ax.set_ylim([0, 2000])
        self.ax.set_xlim([0, self.x_ax_min])
        self.ax.grid(True)

        self.ax.plot([0, self.x_ax_min], [0, 0], 'bo-')

        self.fig.tight_layout()

        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(expand=True, anchor='ne')

        self.__draw_canv__()

    def clear_plot(self):
        """Clear uncertainties from plot to create blank figure (main line is not cleared as this line is updated
        rather than redrawn"""
        num_lines = len(self.stds)
        for i in range(num_lines):
            self.stds[i].pop(0).remove()
        self.stds = []  # Reset list

    def update_plot(self):
        self.ax.lines[0].set_data(self.scan_proc.scan_angles, self.scan_proc.column_densities)

        # Add uncertainty
        self.stds.append(self.ax.plot([self.scan_proc.scan_angles[-1], self.scan_proc.scan_angles[-1]],
                                    [self.doas_worker.column_density['SO2'] - self.doas_worker.std_err,
                                    self.doas_worker.column_density['SO2'] + self.doas_worker.std_err], color='red'))

        # Set x and y limits if we have more than one data point in the array
        if len(self.scan_proc.column_densities) > 1:
            ymax = np.amax(self.scan_proc.column_densities) * 1.15
            ymin = np.amin(self.scan_proc.column_densities)
            if ymin > 0:
                ymin = 0
            else:
                ymin *= 1.1
            if ymax == 0:
                ymax = 2000
            self.ax.set_ylim([ymin, ymax])

            if self.scan_proc.scan_angles[-1] > self.x_ax_min:
                self.ax.set_xlim([0, self.scan_proc.scan_angles[-1] * 1.1])

        self.Q.put(1)

    @property
    def plume_speed(self):
        return self._plume_speed.get()

    @plume_speed.setter
    def plume_speed(self, value):
        self._plume_speed.set(value)

    @property
    def plume_dist(self):
        return self._plume_dist.get()

    @plume_dist.setter
    def plume_dist(self, value):
        self._plume_dist.set(value)

    def update_emission_rate(self):
        """Updates emisison rate label"""
        self.emission_rate.configure(text='{:.2f}'.format(self.scan_proc.SO2_flux))
        self.emission_rate_td.configure(text='{:.1f}'.format(self.scan_proc.flux_tons))

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            if update == 1:
                self.canv.draw()
            else:
                return
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

    def close_widget(self):
        """Closes widget cleanly, by stopping __draw_canv__()"""
        self.Q.put(2)


class TimeSeriesPlot:
    """
    Class to plot column densities retrieved from a scan or traverse sequence
    """
    def __init__(self, root, frame, doas_worker=DOASWorker(), fig_size=(10,3), dpi=96):
        self.root = root

        self.doas_worker = doas_worker
        self.doas_worker.fig_series = self
        self.series = self.doas_worker.series

        self.formatter = DateFormatter('%H:%M:%S')
        self.locator = HourLocator(interval=1)

        self.Q = queue.Queue()

        self.fig_size = fig_size
        self.dpi = dpi

        self.x_ax_min = 20

        self.__setup_gui__(frame)

    def __setup_gui__(self,frame):
        """Controls widget setup"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=5)

        # INPUT FRAME SETUP
        self.frame_inputs = ttk.Frame(self.frame)
        self.frame_inputs.pack(side=tk.LEFT, anchor='nw')

        label = ttk.Label(self.frame_inputs, text='Mean emission rate [kg/s]:').grid(
            row=2, column=0, padx=5, pady=5, sticky='e')
        self.emission_rate_mean = ttk.Label(self.frame_inputs, text='N/A')
        self.emission_rate_mean.grid(row=2, column=1, padx=5, pady=5)

        label = ttk.Label(self.frame_inputs, text='Max. emission rate [kg/s]:').grid(
            row=3, column=0, padx=5, pady=5, sticky='e')
        self.emission_rate_max = ttk.Label(self.frame_inputs, text='N/A')
        self.emission_rate_max.grid(row=3, column=1, padx=5, pady=5)

        # FIGURE SETUP
        self.fig = plt.Figure(figsize=self.fig_size, dpi=self.dpi)
        self.ax = self.fig.subplots(1, 1)
        self.ax.set_ylabel('Emission rate [kg/s]')
        self.ax.set_xlabel('Time')
        self.ax.set_ylim([0, 10])
        start_vals = [date2num(datetime.datetime(2020, 1, 1, 0)), date2num(datetime.datetime(2020, 1, 1, 12))]
        self.ax.set_xlim(start_vals)
        self.ax.grid(True)

        self.ax.plot_date(start_vals, [0, 0], fmt='bo-', xdate=True)
        self.ax.xaxis.set_major_locator(self.locator)
        self.ax.xaxis.set_major_formatter(self.formatter)
        self.ax.xaxis.set_tick_params(rotation=30, labelsize=10)

        self.fig.tight_layout()

        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(expand=True, anchor='ne')

        self.__draw_canv__()

    def update_emission_rate(self):
        """Updates emisison rate label"""
        self.emission_rate_mean.configure(text='{:.2f}'.format(self.series.mean_emission))
        self.emission_rate_max.configure(text='{:.1f}'.format(self.series.max_emission))

    def update_plot(self):
        self.ax.lines[0].set_data(self.series.matplotlib_times, self.series.emission_rates)

        # Clear old error bars before plotting new
        if hasattr(self, 'error_bars'):
            for line in self.error_bars[2]:
                line.remove()
        self.error_bars = self.ax.errorbar(self.series.matplotlib_times, self.series.emission_rates,
                                           yerr=self.series.emission_uncertainties, fmt='none', ecolor='red')

        # Set x and y limits if we have more than one data point in the array
        if len(self.series.emission_rates) > 0:
            ymax = np.amax(self.series.emission_rates)
            ymin = np.amin(self.series.emission_rates)
            if ymin > 0:
                ymin = 0
            else:
                ymin *= 1.2
            if ymax == 0:
                ymax = 1
            elif ymax < 0:
                ymax *= 0.8
            else:
                ymax *= 1.2
            self.ax.set_ylim([ymin, ymax])

            time_delta = datetime.timedelta(minutes=30)
            min = date2num(np.nanmin(self.series.times) - time_delta)
            max = date2num(np.nanmax(self.series.times) + time_delta)
            self.ax.set_xlim([min, max])

        self.Q.put(1)

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            if update == 1:
                self.canv.draw()
            else:
                return
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

    def close_widget(self):
        """Closes widget cleanly, by stopping __draw_canv__()"""
        self.Q.put(2)