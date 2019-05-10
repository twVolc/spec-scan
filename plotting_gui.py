from matplotlib import pyplot as plt

from gui_subs import SettingsGUI

class SpectraPlot:
    """
    Generates a widget containing 3 subplots of spectra -> dark, clear (Fraunhofer), in-plume
    """

    def __init__(self, frame, doas_worker):
        self.setts = SettingsGUI()
        self.doas_worker = doas_worker

    def __setup_gui__(self):
        pass
