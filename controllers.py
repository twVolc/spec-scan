# IMPORTING SEABREEZE WITH PYUSB BACKEND FOR SPECTROMETER CONTROL
import seabreeze.spectrometers as sb


class SpecCtrl:
    """Class to control spectrometer acquisition from USB2000+/Flame spectrometers"""

    def __init__(self, int_time=100000):
        # Discover spectrometer devices
        self.devices = sb.list_devices()
        self.spec = sb.Spectrometer(self.devices[0])

        # Set integration time (ALL IN MICROSECONDS)
        self._int_limit_lower = 1000        # Lower integration time limit
        self._int_limit_upper = 6000000     # Upper integration time limit
        self._int_time = None               # Integration time attribute
        self.int_time = int_time

    @property
    def int_time(self):
        return self._int_time

    @int_time.setter
    def int_time(self, int_time):
        """Set integration time"""
        # Check requested integration time is acceptable
        if int_time < self._int_limit_lower:
            raise ValueError('Integration time below %i us is not possible' % self._int_limit_lower)
        elif int_time > self._int_limit_upper:
            raise ValueError('Integration time above %i us is not possible' % self._int_limit_upper)

        self._int_time = int_time
        self.spec.integration_time_micros(int_time)

    def get_spec(self):
        """Acquire spectrum from spectrometer"""
        return self.spec.intensities()

    def get_wavelengths(self):
        """Returns wavelengths"""
        return self.spec.wavelengths()


class AcqController:
    """Handles overall control of spectral acquisition, including motor movement"""
    def __init__(self):
        self.scan_deg = None    # Degrees of scan







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
