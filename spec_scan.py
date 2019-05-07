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
        self.mainFrame = ttk.Frame(self.tabs)
        self.calibFrame = ttk.Frame(self.tabs)
        self.hyperspecFrame = ttk.Frame(self.tabs)
        self.tabs.add(self.mainFrame, text='Main')
        self.tabs.add(self.calibFrame, text='Calibration')
        self.tabs.add(self.hyperspecFrame, text='Hyperspectral')
        self.tabs.pack(fill="both", expand=1)

        # ==============================================================================================================
        # MAIN TAB SETUP
        # ==============================================================================================================
        self.messages = MessagesGUI(self.mainFrame)
        self.messages.frame.pack(side='right', fill=tk.Y, expand=1, anchor='e')

        # ==============================================================================================================
        # Calibration work - reference spectrum etc
        # ==============================================================================================================
        # -------------------------
        # Reference Spectrum Setup
        # -------------------------
        self.refFrame = tk.LabelFrame(self.calibFrame, text='Reference Spectrum', font=self.mainFontBold,
                                      relief=tk.RAISED, borderwidth=2, bg=self.bgColour)
        self.refFrame.grid(row=0, column=0, padx=self.pdx, pady=self.pdy)

        self.loadRefFrame = tk.Frame(self.refFrame)
        self.loadRefFrame.grid(row=0, column=0, sticky='w')
        self.labelRef = tk.Label(self.loadRefFrame, text='Filename:', font=self.mainFontBold)
        self.labelRef.grid(row=0, column=0, padx=self.pdx, pady=self.pdy)
        self.nameRef = tk.Label(self.loadRefFrame, text='None Selected', font=self.mainFont)
        self.nameRef.grid(row=0, column=1, padx=self.pdx, pady=self.pdy)
        self.selectRef = tk.Button(self.loadRefFrame, text='Load Spectrum', command=self.choose_ref_spec)
        self.selectRef.grid(row=0, column=2, padx=self.pdx, pady=self.pdy)

        # Plot reference spectrum
        self.FigRef = plt.Figure(figsize=self.ref_fig_size, dpi=self.dpi)
        self.AxRef = self.FigRef.add_subplot(111)

        self.AxRef.set_title('Reference spectrum: None')
        self.AxRef.set_ylabel(r'Absorption Cross Section [cm$^2$/molecule]')
        self.AxRef.set_xlabel('Wavelength [nm]')
        self.AxRef.tick_params(axis='both', direction='in', top='on', right='on')

        self.refCanvas = FigureCanvasTkAgg(self.FigRef, master=self.refFrame)
        self.refCanvas.draw()
        self.refCanvas.get_tk_widget().grid(row=1, column=0)
        # --------------------------------------------------------------------------------------------------------------

        # -----------------------
        # Calibration image setup
        # -----------------------
        self.calFrame = tk.LabelFrame(self.calibFrame, text='PiSpec Calibration Spectrum', font=self.mainFontBold,
                                      relief=tk.RAISED, borderwidth=2, bg=self.bgColour)
        self.calFrame.grid(row=0, column=1, padx=self.pdx, pady=self.pdy)

        self.loadCalFrame = tk.Frame(self.calFrame)
        self.loadCalFrame.grid(row=0, column=0, sticky='w')
        self.labelCal = tk.Label(self.loadCalFrame, text='Filename:', font=self.mainFontBold)
        self.labelCal.grid(row=0, column=0, padx=self.pdx, pady=self.pdy)
        self.nameCal = tk.Label(self.loadCalFrame, text='None Selected', font=self.mainFont)
        self.nameCal.grid(row=0, column=1, padx=self.pdx, pady=self.pdy)
        self.selectCal = tk.Button(self.loadCalFrame, text='Load Image', command=self.choose_cal_img)
        self.selectCal.grid(row=0, column=2, padx=self.pdx, pady=self.pdy)

        # Plot reference spectrum
        fig_height_ratios = [5, 2]
        self.AxCal = [0, 0]
        self.FigCal = plt.Figure(figsize=self.cal_fig_size, dpi=self.dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=fig_height_ratios)
        # gs.update(wspace= 0.000000001, hspace=0.00000000001)
        self.AxCal[0] = self.FigCal.add_subplot(gs[0])
        self.AxCal[0].set_aspect(1)
        self.AxCal[1] = self.FigCal.add_subplot(gs[1])
        self.AxCal[1].set_aspect((fig_height_ratios[1]/fig_height_ratios[0]) * (self.imgSizeY/self.maxDN))
        self.FigCal.subplots_adjust(hspace=0.1)
        self.AxCal[0].tick_params(axis='both', direction='in', top='on', right='on')
        self.AxCal[1].tick_params(axis='both', direction='in', top='on', right='on')

        self.AxCal[0].set_title('Calibration spectrum: None')
        self.AxCal[0].set_ylabel('Pixel')
        self.AxCal[0].set_ylim([self.imgSizeY-1, 0])
        self.AxCal[0].set_xlim([0, self.imgSizeX-1])
        self.AxCal[0].set_xticklabels([])


        self.AxCal[1].set_ylabel('DN')
        self.AxCal[1].set_xlabel('Pixel')
        self.AxCal[1].set_ylim([0, self.maxDN])
        self.AxCal[1].set_xlim([0, self.imgSizeX-1])

        self.calCanvas = FigureCanvasTkAgg(self.FigCal, master=self.calFrame)
        self.calCanvas.draw()
        self.calCanvas.get_tk_widget().grid(row=1, column=0)



    def choose_ref_spec(self):
        """Load reference spectrum"""
        self.ref_spec_path = filedialog.askopenfilename(initialdir=self.ref_spec_init_dir,
                                                        title='Select reference spectrum',
                                                        filetypes=(("Text files", "*.txt"), ("All files", "*.*")))
        self.ref_spec_file = self.ref_spec_path.split('/')[-1]
        if not self.ref_spec_path:
            return
        if len(self.ref_spec_path) > 53:
            self.nameRef.configure(text='...' + self.ref_spec_path[-50:])
        else:
            self.nameRef.configure(text=self.ref_spec_path)
        self.DOAS.load_reference_spectrum(self.ref_spec_path)

        self.plot_ref_spec()

    def plot_ref_spec(self):
        """Plot up reference spectrum"""
        if hasattr(self, 'ref_plot'):  # If we have already plotted a spectra we can just update
            self.ref_plot.set_data(self.DOAS.ref_spec[:, 0], self.DOAS.ref_spec[:, 1])
        else:
            self.ref_plot, = self.AxRef.plot(self.DOAS.ref_spec[:, 0], self.DOAS.ref_spec[:, 1])
        self.AxRef.set_xlim([self.DOAS.ref_spec[0, 0], self.DOAS.ref_spec[-1, 0]])
        self.AxRef.set_ylim([0, np.amax(self.DOAS.ref_spec[:, 1])])
        self.AxRef.set_title('Reference Spectrum: %s' % self.ref_spec_file)
        self.refCanvas.draw()

    def choose_cal_img(self):
        """Load reference spectrum"""
        self.cal_spec_path = filedialog.askopenfilename(initialdir=self.local_dir,
                                                        title='Select calibration image',
                                                        filetypes=(("png files", "*.png"), ("All files", "*.*")))
        self.cal_spec_file = self.cal_spec_path.split('/')[-1]
        if not self.cal_spec_path:
            return
        if len(self.cal_spec_path) > 53:
            self.nameCal.configure(text='...' + self.cal_spec_path[-50:])
        else:
            self.nameCal.configure(text=self.cal_spec_path)
        self.DOAS.load_calibration_spectrum(self.cal_spec_path)


        self.cal_plot_ctrl()

    def cal_plot_ctrl(self):
        """Control the plotting of both calibration figures -> draw canvas"""
        self.plot_cal_img()
        self.plot_cal_spec()
        self.calCanvas.draw()

    def plot_cal_img(self):
        """Display calibration image"""
        if hasattr(self, 'cal_img_plot'):  # Plot calibation image
            self.cal_img_plot.set_data(self.DOAS.cal_img)
        else:
            self.cal_img_plot = self.AxCal[0].imshow(self.DOAS.cal_img, cmap=cm.Greys_r, interpolation='none', vmin=0,
                                                    vmax=np.amax(self.maxDN))
        self.AxCal[0].set_title('Calibration spectrum: %s' % self.cal_spec_file)

    def plot_cal_spec(self):
        """Display calibration spectrum"""
        if hasattr(self, 'cal_spec_plot'):
            self.cal_spec_plot.set_data(self.DOAS.pixel_vector, self.DOAS.cal_spec)
        else:
            self.cal_spec_plot, = self.AxCal[1].plot(self.DOAS.pixel_vector, self.DOAS.cal_spec)


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