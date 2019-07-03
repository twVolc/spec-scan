# GUI for DOAS processing

# GUI libraries
import tkinter as tk
from tkinter import messagebox
import tkinter.ttk as ttk
import tkinter.font as tkFont
from tkinter import filedialog
from PIL import Image, ImageTk

# Matplotlib dependencies
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets     # For rectangle crop
import matplotlib.patches as patches     # For rectangle crop
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.dates as mdates        # For plotting time on x-axis of flux time series
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler  # implement the default mpl key bindings
from matplotlib.figure import Figure
# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})  # Can make figure labels visible but also does strange things

import numpy as np
import sys

# In-house modules
from doas_routine import DOASWorker, ScanProcess
from config_parser import config_parser
from gui_subs import *
from acquisition_gui import AcquisitionFrame
from postprocess_gui import PostProcess
from plotting_gui import SpectraPlot, DOASPlot, CDPlot
from calibration_gui import CalPlot, RefPlot

class PySpec(ttk.Frame):
    '''PySpec GUI'''
    def __init__(self, parent, x_size, y_size):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        parent.title('PySpec')
        self.parent.protocol('WM_DELETE_WINDOW', self.exit_app)

        self.version = '2.0, September 17 2017.'  # PyCam Version

        self.DOAS = DOASWorker(1)
        self.scan_proc = ScanProcess()

        self.config = config_parser()

        # ==============================================================================================================
        # GUI SETUP
        # ==============================================================================================================
        matplotlib.rcParams.update({'font.size': self.config['mtplt_font_size']})
        self.bgColour = "#ccc"
        style = ttk.Style()
        style.configure("TFrame", background=self.bgColour)

        self.pdx = 2
        self.pdy = 2

        ratio_x = x_size / 1920  # Find ratio of screen resolution to my screen
        ratio_y = y_size / 1080
        self.ratio = min([ratio_x, ratio_y])
        self.font_size = int(np.round(10 * self.ratio))  # Calculate a suitable fontsize
        self.mainFont = tkFont.Font(family='arial', size=self.font_size)
        self.mainFontBold = tkFont.Font(family='arial', size=self.font_size, weight='bold')

        self.tabs = ttk.Notebook(self.parent)
        self.mainFrame = ttk.Frame(self.tabs, borderwidth=2)
        self.calibFrame = ttk.Frame(self.tabs, borderwidth=2)
        self.tabs.add(self.mainFrame, text='Main')
        self.tabs.add(self.calibFrame, text='Calibration')
        self.tabs.pack(fill="both", expand=1)

        # ==============================================================================================================
        # MAIN TAB SETUP
        # ==============================================================================================================
        # Scrolling canvas setup for main frame
        self.main_canvas = tk.Canvas(self.mainFrame, borderwidth=0, background=self.bgColour)
        self.main_canvas_scroll = ScrollWindow(self.mainFrame, self.main_canvas)
        self.frame_1 = ttk.Frame(self.main_canvas_scroll.frame, borderwidth=2)
        self.frame_1.pack(expand=True, fill=tk.BOTH, anchor='nw')

        # Messages frame
        self.messages = MessagesGUI(self.frame_1)
        self.messages.frame.pack(side='right', fill=tk.Y, anchor='e', expand=False)

        # Plot frame
        self.plot_frame = ttk.Frame(self.frame_1)
        self.plot_frame.pack(side='right', expand=1, anchor='e')

        # Spectra + DOAS plots
        self.doas_frame = DOASPlot(self.parent, self.plot_frame, self.DOAS,
                                   figsize=self.config['doas_fig_size'], dpi=self.config['dpi'])
        self.spec_frame = SpectraPlot(self.parent, self.plot_frame, self.DOAS, self.doas_frame,
                                      figsize=self.config['spec_fig_size'], dpi=self.config['dpi'])
        self.spec_frame.frame.pack(side='top', expand=1, anchor='n', fill=tk.X)
        self.doas_frame.frame.pack(side='top', expand=1, anchor='n', fill=tk.X)

        self.cd_plot = CDPlot(self.parent, self.plot_frame, scan_proc=self.scan_proc,
                              fig_size=self.config['scan_fig_size'], dpi=self.config['dpi'])
        self.cd_plot.frame.pack(side='top', fill=tk.BOTH, expand=1, anchor='nw')

        # Left side of main tab GUI
        self.left_main_frame = ttk.Frame(self.frame_1)
        self.left_main_frame.pack(side='left', expand=1, fill=tk.BOTH)

        # Acquisition frame
        self.acq_frame = AcquisitionFrame(self.left_main_frame, self.DOAS, self.scan_proc,
                                          self.spec_frame, self.doas_frame, self.cd_plot, self.config['arduino_COM'],
                                          save_path=self.config['init_dir'])
        self.acq_frame.frame.pack(side='top', expand=False, anchor='nw', fill=tk.X)
        self.doas_frame.acq_obj = self.acq_frame

        self.post_process_frame = PostProcess(self.left_main_frame, self.DOAS, self.scan_proc,
                                              self.spec_frame, self.doas_frame, self.cd_plot)
        self.post_process_frame.frame.pack(side='top', expand=True, anchor='nw', fill=tk.X)

        # ==============================================================================================================
        # Calibration work - reference spectrum etc
        # ==============================================================================================================
        self.cal_canvas = tk.Canvas(self.calibFrame, borderwidth=0, background=self.bgColour)
        self.cal_canvas_scroll = ScrollWindow(self.calibFrame, self.cal_canvas)
        self.frame_2 = ttk.Frame(self.cal_canvas_scroll.frame, borderwidth=2)
        self.frame_2.pack(expand=True, fill=tk.BOTH, anchor='nw')

        self.ILS_frame = CalPlot(self.frame_2, self.acq_frame.spec_ctrl, self.DOAS, config=self.config)
        self.ILS_frame.frame.pack(side=tk.RIGHT, fill=tk.Y, expand=1, anchor='e')

        self.ref_frame = RefPlot(self.frame_2, self.DOAS, self.config['ref_spec_dir'], self.config['ref_fig_size'],
                                 self.config['dpi'], self.config['ref_spec_SO2'])
        self.ref_frame.frame.pack(side=tk.LEFT, anchor='n', padx=5, pady=5)


    def exit_app(self):
        """Exit GUI"""
        if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
            self.parent.destroy()
            sys.exit()


def run_GUI():
    padx = 0
    pady = 50
    root = tk.Tk()
    root.geometry('{0}x{1}+0+0'.format(root.winfo_screenwidth() - padx, root.winfo_screenheight() - pady))
    x_size = root.winfo_screenwidth()  # Get screen width
    y_size = root.winfo_screenheight()  # Get screen height
    myGUI = PySpec(root, x_size, y_size)
    root.mainloop()

if __name__ == '__main__':
    run_GUI()