# IMPORTING SEABREEZE WITH PYUSB BACKEND FOR SPECTROMETER CONTROL
import seabreeze.spectrometers as sb
import numpy as np


class SpecCtrl:
    """
    Class to control spectrometer acquisition from USB2000+/Flame spectrometers
    """
    def __init__(self, int_time=100, coadd=1):
        # Discover spectrometer devices
        self.devices = None     # List of detected spectrometers
        self.spec = None        # Holds spectrometer for interfacing via seabreeze
        self.wavelengths = None # Array of wavelengths
        self.find_device()

        # Set integration time (ALL IN MICROSECONDS)
        self._int_limit_lower = 1000        # Lower integration time limit
        self._int_limit_upper = 20000000    # Upper integration time limit
        self._int_time = None               # Integration time attribute
        self.int_time = int_time

        self.min_coadd = 1
        self.max_coadd = 100
        self._coadd = None      # Controls coadding of spectra
        self.coadd = coadd


    def find_device(self):
        """Function to search for devices"""
        try:
            self.devices = sb.list_devices()
            self.spec = sb.Spectrometer(self.devices[0])
            self.spec.trigger_mode(0)

            # If we have a spectrometer we then retrieve its wavelength calibration and store it as an attribute
            self.get_wavelengths()

        except IndexError:
            self.devices = None
            self.spec = None
            raise SpectrometerConnectionError('No spectrometer found')


    @property
    def int_time(self):
        return self._int_time / 1000  # Return time in milliseconds

    @int_time.setter
    def int_time(self, int_time):
        """Set integration time"""
        int_time *= int(1000)       # Adjust to work in microseconds (class takes time in milliseconds)

        # Check requested integration time is acceptable
        if int_time < self._int_limit_lower:
            raise ValueError('Integration time below %i us is not possible' % self._int_limit_lower)
        elif int_time > self._int_limit_upper:
            raise ValueError('Integration time above %i us is not possible' % self._int_limit_upper)

        self._int_time = int_time
        self.spec.integration_time_micros(int_time)

    @property
    def coadd(self):
        return self._coadd

    @coadd.setter
    def coadd(self, coadd):
        """Set coadding property"""
        if coadd < self.min_coadd:
            coadd = self.min_coadd
        elif coadd > self.max_coadd:
            coadd = self.max_coadd
        coadd = int(coadd)
        self._coadd = coadd

    def get_spec(self):
        """Acquire spectrum from spectrometer"""
        # First spectrum is discarded as it may have been partially acquired prior to the acquisition request
        # Trigger_mode is continuous
        _dummy = self.spec.intensities()

        # Set array for coadding spectra
        coadded_spectrum = np.zeros(len(self.wavelengths))

        # Loop through number of coadds
        for i in range(self.coadd):
            coadded_spectrum += self.spec.intensities()

        # Correct for number of coadds to result in a spectrum with correct digital numbers for bit-depth of device
        coadded_spectrum /= self.coadd
        return coadded_spectrum

    def get_spec_now(self):
        """Immediately acquire spectrum from spectrometer - does not discard first spectrum (probably never used)"""
        return self.spec.intensities()

    def get_wavelengths(self):
        """Returns wavelengths"""
        self.wavelengths = self.spec.wavelengths()


class SpectrometerConnectionError(Exception):
    """
    Error raised if no spectrometer is detected
    """
    pass


class ScanProperties:
    """
    Handles overall control of spectral acquisition, including motor movement
    """
    def __init__(self):
        self.scan_incr = 1.8  # Stepper motors scan increment
        self.scan_range_full = 100  # Range of motion in scanner
        self.return_steps = int(
            np.ceil(self.scan_range_full / self.scan_incr))  # Max num of steps to reset motor at start
        self.scan_fwd = b'\x00'  # Set here whether 0 or 1 steps the scanner in the forward direction
        self.scan_back = b'\x01'

        self.dark_steps = 50    # Number of steps to get to dark position


if __name__=="__main__":

    from matplotlib import pyplot as plt
    import matplotlib.animation as animation

    spec = SpecCtrl()
    spec.int_time = 100000

    wavelength = spec.get_wavelengths()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    def animate(i):
        spectrum = spec.get_spec()
        ax.clear()
        ax.plot(wavelength, spectrum)

    ani = animation.FuncAnimation(fig, animate, interval=100)
    plt.show()
