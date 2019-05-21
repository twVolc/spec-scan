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
from doas_routine import DOASWorker
from config_parser import config_parser
from gui_subs import *
from acquisition_gui import AcquisitionFrame
from plotting_gui import SpectraPlot, DOASPlot
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

        self.local_dir, self.ref_spec_init_dir, self.ref_fig_size, self.cal_fig_size, \
        self.dpi, mtplt_font_size, self.imgSizeX, self.imgSizeY = config_parser()

        matplotlib.rcParams.update({'font.size': mtplt_font_size})

        self.maxDN = (2 ** 10) - 1  # Maximum digital number of images/spectra

        # ==============================================================================================================
        # GUI SETUP
        # ==============================================================================================================
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
        # Messages frame
        self.messages = MessagesGUI(self.mainFrame)
        self.messages.frame.pack(side='right', fill=tk.Y, expand=1, anchor='e')

        # Plot frame
        self.plot_frame = ttk.Frame(self.mainFrame)
        self.plot_frame.pack(side='right', expand=1, anchor='e')

        # Spectra plots
        self.spec_frame = SpectraPlot(self.plot_frame, self.DOAS)
        self.spec_frame.frame.pack(side='top', expand=1, anchor='n')

        # DOAS plot
        self.doas_frame = DOASPlot(self.plot_frame, self.DOAS)
        self.doas_frame.frame.pack(side='top', expand=1, anchor='n')

        # Acquisition frame
        self.acq_frame = AcquisitionFrame(self.mainFrame, self.DOAS, self.spec_frame)
        self.acq_frame.frame.pack(side='left', expand=1, anchor='nw')

        # ==============================================================================================================
        # Calibration work - reference spectrum etc
        # ==============================================================================================================
        self.ILS_frame = CalPlot(self.calibFrame, self.acq_frame.spec_ctrl, self.DOAS)
        self.ILS_frame.frame.pack(side=tk.RIGHT, fill=tk.Y, expand=1, anchor='e')

        self.ref_frame = RefPlot(self.calibFrame, self.DOAS, self.ref_spec_init_dir, self.ref_fig_size, self.dpi)
        self.ref_frame.frame.pack(side=tk.LEFT)

        

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