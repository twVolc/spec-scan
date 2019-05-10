import tkinter as tk
import tkinter.ttk as ttk

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import queue

from gui_subs import SettingsGUI
from doas_routine import DOASWorker

class SpectraPlot:
    """
    Generates a widget containing 3 subplots of spectra -> dark, clear (Fraunhofer), in-plume
    """

    def __init__(self, frame, doas_worker=DOASWorker()):
        self.setts = SettingsGUI()
        self.doas_worker = doas_worker

        self.max_DN = 2**16 - 1  # Maximum DN for spectrometer

        # Could use threads and queues to update plots, or just us simple functions which Acquisition frame calls
        self.dark_Q = queue.Queue()
        self.clear_Q = queue.Queue()
        self.plume_Q = queue.Queue()

        self.__setup_gui__(frame)

    def __setup_gui__(self, frame):
        """Organise widget"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=4)

        # -----------------------------------------------
        # WIDGET SETUP
        # ------------------------------------------------
        self.frame2 = ttk.Frame(self.frame)
        self.frame2.pack(side=tk.TOP, fill=tk.X, expand=1)

        label = tk.Label(self.frame2, text='Fit wavelength (min.):').pack(side=tk.LEFT)
        # label.grid(row=0, column=0)
        self.fit_wind_start = tk.DoubleVar()
        self.fit_wind_start.set(self.doas_worker.start_fit_wave)
        self.fit_wind_box_start = tk.Spinbox(self.frame2, from_=0, to=400, increment=0.1, width=5,
                                             textvariable=self.fit_wind_start, command=self.update_fit_wind_start)
        # self.fit_wind_box_start.grid(row=0, column=1)
        self.fit_wind_box_start.pack(side=tk.LEFT)

        label2 = tk.Label(self.frame2, text='Fit wavelength (max.):').pack(side=tk.LEFT)
        # label2.grid(row=0, column=2)
        self.fit_wind_end = tk.DoubleVar()
        self.fit_wind_end.set(self.doas_worker.end_fit_wave)
        self.fit_wind_box_end = tk.Spinbox(self.frame2, from_=1, to=400, increment=0.1, width=5,
                                           textvariable=self.fit_wind_end, command=self.update_fit_wind_end)
        # self.fit_wind_box_end.grid(row=0, column=3)
        self.fit_wind_box_end.pack(side=tk.LEFT)

        # ------------------------------------------------
        # FIGURE SETUP
        # ------------------------------------------------
        plt.style.use('dark_background')
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

        self.min_line = self.ax.plot([self.doas_worker.start_fit_wave,self.doas_worker.start_fit_wave],
                                     [0, self.max_DN], 'y')
        self.max_line = self.ax.plot([self.doas_worker.end_fit_wave,self.doas_worker.end_fit_wave],
                                     [0, self.max_DN], 'y')

        width = self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave
        self.fit_wind = Rectangle((self.doas_worker.start_fit_wave,0), width=width, height=self.max_DN,
                                  fill=True, facecolor='y', alpha=0.2)
        self.fit_wind_patch = self.ax.add_patch(self.fit_wind)

        self.canv = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canv.draw()
        self.canv.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def update_dark(self):
        """Update dark plot with new data"""
        self.ax.lines[0].set_data(self.doas_worker.wavelengths, self.doas_worker.dark_spec)
        self.ax.set_xlim([self.doas_worker.wavelengths[0],self.doas_worker.wavelengths[-1]])
        self.canv.draw()

    def update_clear(self):
        """Update clear plot with new data"""
        self.ax.lines[1].set_data(self.doas_worker.wavelengths, self.doas_worker.clear_spec_raw)
        self.ax.set_xlim([self.doas_worker.wavelengths[0], self.doas_worker.wavelengths[-1]])
        self.canv.draw()

    def update_plume(self):
        """Update clear plot with new data"""
        self.ax.lines[2].set_data(self.doas_worker.wavelengths, self.doas_worker.plume_spec_raw)
        self.ax.set_xlim([self.doas_worker.wavelengths[0], self.doas_worker.wavelengths[-1]])
        self.canv.draw()

    def update_fit_wind_start(self):
        """updates fit window on plot"""
        # Update DOASWorker with new fit window
        self.doas_worker.start_fit_wave = self.fit_wind_start.get()

        # Update plot
        self.min_line[0].set_data([self.doas_worker.start_fit_wave, self.doas_worker.start_fit_wave], [0, self.max_DN])
        self.fit_wind.set_x(self.doas_worker.start_fit_wave)
        self.fit_wind.set_width(self.doas_worker.end_fit_wave-self.doas_worker.start_fit_wave)
        self.canv.draw()

    def update_fit_wind_end(self):
        """updates fit window on plot"""
        # Update DOASWorker with new fit window
        self.doas_worker.end_fit_wave = self.fit_wind_end.get()

        # Update plot
        self.max_line[0].set_data([self.doas_worker.end_fit_wave, self.doas_worker.end_fit_wave], [0, self.max_DN])
        self.fit_wind.set_width(self.doas_worker.end_fit_wave - self.doas_worker.start_fit_wave)
        self.canv.draw()


class DOASPlot:
    """
    Generates a widget containing the DOAS fit plot
    """
    def __init__(self, frame, doas_worker=DOASWorker()):
        self.setts = SettingsGUI()
        self.doas_worker = doas_worker

    def __setup_gui__(self, frame):
        """Organise widget"""
        self.frame = ttk.Frame(frame, relief=tk.RAISED, borderwidth=4)
