# -*- coding: utf-8 -*-

"""
Script to be run for carrying out automatic DOAS scans
"""

# Import standard library packages
import sys
import datetime
import os

# Enter absolute path to this package on the filesystem
path_to_package = '../'
sys.path.append(path_to_package)

# Import custom packages
from automation.controllers import Spectrometer, Scanner
from automation.utils import format_time, save_spectrum

# -------------------------------------------
# Define path to where spectra will be saved
# I would recommend mounting an external hard drive to the Raspberry Pi and pointing this save path there
# This will protect from data loss if external data transmission breaks down
spectra_dir = './Spectra/'
# -------------------------------------------

# Check that Spectra directory exists, if not - make it
if os.path.exists(spectra_dir):
    os.mkdir(spectra_dir)

# Insantiate our spectrometer and scanner objects
spectrometer = Spectrometer()
scanner = Scanner()

# Reset scanner to dark position
scanner.full_return()

# Continuous loop through scanning
saturation_list = []    # List containing saturation levels within spectra
while True:

    # Check saturation levels - we first check for over saturation, as this is critical to spectrum integrity
    # We then check for undersaturation which is not vital but will cause noise issues
    if -1 in saturation_list:
        try:
            spectrometer.int_time_idx -= 1  # Updating the int_time_idx automatically updates int_time
        except IndexError:
            pass
    elif 1 in saturation_list:
        try:
            spectrometer.int_time_idx += 1
        except IndexError:
            pass

    # Reset saturation list
    saturation_list = []

    # Get time in formatted string
    time_str = format_time(datetime.datetime.now(), spectrometer.file_datestr)

    # Make a new directory for this scan
    scan_dir = spectra_dir + time_str + '/'
    os.mkdir(scan_dir)

    # Create a lock indicating that this scan is not finished yet
    lock = scan_dir + 'scan.lock'
    open(lock, 'w').close()

    # Acquire dark spectrum
    spectrometer.get_spec()

    # Generate filename
    filename = spectrometer.generate_filename(time_str, spectrometer.file_spec_type['dark'])

    # Save spectrum
    save_spectrum(spectrometer.wavelengths, spectrometer.spectrum, scan_dir + filename)

    # Step to first clear sky spectrum
    scanner.dark_to_clear()

    # ==================================================================================================================
    # Main loop for acquiring scan spectra
    for i in range(Scanner.scan_properties.total_scan_steps):

        # Get time in formatted string
        time_str = format_time(datetime.datetime.now(), spectrometer.file_datestr)

        # Acquire spectrum
        spectrometer.get_spec()

        # Generate filename
        filename = spectrometer.generate_filename(time_str, spectrometer.file_spec_type['meas'])

        # Save spectrum
        save_spectrum(spectrometer.wavelengths, spectrometer.spectrum, scan_dir + filename)

        # Step motor forwards
        scanner.step_forward()

        # Check saturation of spectrum - append result to saturation list
        saturation_list.append(spectrometer.check_saturation())

    # Remove lock, indicating that this scan is now ready for transmission
    os.remove(lock)

    # Step back to the start of clear sky and then to dark position, to start new scan
    scanner.scan_return()
    scanner.clear_to_dark()


