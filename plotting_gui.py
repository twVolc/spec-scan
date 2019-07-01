import tkinter as tk
import tkinter.ttk as ttk

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np
import queue

from gui_subs import SettingsGUI
from doas_routine import DOASWorker, ScanProcess
from acquisition_gui import AcquisitionFrame

plt.style.use('dark_background')

refresh_rate = 100

class SpectraPlot:
    """
    Generates a widget containing 3 subplots of spectra -> dark, clear (Fraunhofer), in-plume
    """

    def __init__(self, root, frame, doas_worker=DOASWorker(), doas_plot=None):
        self.root = root
        self.setts = SettingsGUI()
        self.doas_worker = doas_worker
        self.doas_plot = doas_plot

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
        self.fig = plt.Figure(figsize=(10, 3), dpi=100)

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

        self.__draw_canv__()

    def update_dark(self):
        """Update dark plot with new data"""
        self.ax.lines[0].set_data(self.doas_worker.wavelengths, self.doas_worker.dark_spec)
        self.ax.set_xlim([self.doas_worker.wavelengths[0],self.doas_worker.wavelengths[-1]])
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
            self.canv.draw()
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

class DOASPlot:
    """
    Generates a widget containing the DOAS fit plot
    """
    def __init__(self, root, frame, doas_worker=DOASWorker()):
        self.root = root

        self.setts = SettingsGUI()
        self.doas_worker = doas_worker

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

        # ------------------------------------------------
        # FIGURE SETUP
        # ------------------------------------------------
        self.fig = plt.Figure(figsize=(10, 3), dpi=100)

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
        self.ax.set_title('Column density [ppm.m]: N/A')
        self.fig.tight_layout()

        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.__draw_canv__()

    def update_shift(self):
        """Updates DOASWorker shift value for aligning spectra"""
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
        # Update plot lines with new data
        self.ax.lines[0].set_data(self.doas_worker.wavelengths_cut, self.doas_worker.abs_spec_cut)
        self.ax.lines[1].set_data(self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit['SO2'])

        # Set axis limits
        self.ax.set_xlim([self.doas_worker.wavelengths_cut[0], self.doas_worker.wavelengths_cut[-1]])
        ylims = np.amax(np.absolute(self.doas_worker.abs_spec_cut))
        ylims *= 1.15
        if ylims == 0:
            ylims = 0.05
        self.ax.set_ylim([-ylims, ylims])
        self.ax.set_title('Column density [ppm.m]: {}'.format(self.doas_worker.column_amount))

        # Draw updates
        self.Q.put(1)

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            self.canv.draw()
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)

    def save_spectra(self):
        """Saves processed DOAS spectra"""
        if isinstance(self.acq_obj, AcquisitionFrame):
            print('Saving...')
            self.acq_obj.save_processed_spec()



class CDPlot:
    """
    Class to plot column densities retrieved from a scan or traverse sequence
    """
    def __init__(self, root, frame, fig_size=(5,5), dpi=100, scan_proc=ScanProcess()):
        self.root = root

        self.Q = queue.Queue()

        self.fig_size=fig_size
        self.dpi = dpi

        self.x_ax_min = 20

        self.scan_proc = scan_proc

        self.__setup_gui__(frame)

    def __setup_gui__(self,frame):
        """Controls widget setup"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=5)

        # INPUT FRAME SETUP
        self.frame_inputs = ttk.Frame(self.frame)
        self.frame_inputs.pack(side=tk.LEFT, anchor='nw')

        label_dist = ttk.Label(self.frame_inputs, text='Plume distance [m]:')
        self.plume_dist = tk.IntVar()
        self.plume_dist.set(1000)
        self.plume_dist_box = ttk.Entry(self.frame_inputs, width=5, textvariable=self.plume_dist)

        label_dist.grid(row=0, column=0, padx=5, pady=5, sticky='e')
        self.plume_dist_box.grid(row=0, column=1, padx=5, pady=5)

        label_speed = ttk.Label(self.frame_inputs, text='Plume speed [m/s]:')
        self.plume_speed = tk.DoubleVar()
        self.plume_speed.set(5.0)
        self.plume_speed_box = ttk.Entry(self.frame_inputs, width=5, textvariable=self.plume_speed)

        label_speed.grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.plume_speed_box.grid(row=1, column=1, padx=5, pady=5)

        label = ttk.Label(self.frame_inputs, text='Scan emission rate [kg/s]:').grid(
            row=2, column=0, padx=5, pady=5, sticky='e')
        self.emission_rate = ttk.Label(self.frame_inputs, text='N/A')
        self.emission_rate.grid(row=2, column=1, padx=5, pady=5)

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

    def update_plot(self):
        self.ax.lines[0].set_data(self.scan_proc.scan_angles, self.scan_proc.column_densities)

        # Set x and y limits if we have more than one data point in the array
        if len(self.scan_proc.column_densities) > 1:
            ymax = np.amax(self.scan_proc.column_densities) * 1.15
            ymin = np.amin(self.scan_proc.column_densities)
            if ymin > 0:
                ymin = 0
            else:
                ymin *= 1.1
            if ymax == 0:
                ymax == 2000
            self.ax.set_ylim([ymin, ymax])

            if self.scan_proc.scan_angles[-1] > self.x_ax_min:
                self.ax.set_xlim([0, self.scan_proc.scan_angles[-1] * 1.1])

        self.Q.put(1)

    def update_emission_rate(self):
        """Updates emisison rate label"""
        self.emission_rate.configure(text='{:.2f}'.format(self.scan_proc.SO2_flux))

    def __draw_canv__(self):
        """Draws canvas periodically"""
        try:
            update = self.Q.get(block=False)
            self.canv.draw()
        except queue.Empty:
            pass
        self.root.after(refresh_rate, self.__draw_canv__)
