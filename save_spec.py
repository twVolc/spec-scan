import numpy as np

from doas_routine import *
from acquisition_gui import *
from controllers import SpecCtrl


class SaveSpectra:
    """
    Handles all saving of spectra from spectrometer acquisitions and DOAS retrievals
    """
    def __init__(self, doas=DOASWorker(), scan_proc=ScanProcess(), spec_ctrl=SpecCtrl()):
        self.doas_worker = doas
        self.spec_ctrl = spec_ctrl
        self.scan_proc = scan_proc

        self.scan_filename = 'Scan_data.txt'


    def save_dark(self, filename):
        """Save dark spectrum"""
        # Sanity checks
        if self.doas_worker.wavelengths is None or self.doas_worker.dark_spec is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.doas_worker.wavelengths) != len(self.doas_worker.dark_spec):
            raise ValueError('Arrays are not the same length. Cannot save.')

        # Save file
        np.savetxt(filename, np.transpose([self.doas_worker.wavelengths, self.doas_worker.dark_spec]),
                   header='Dark spectrum\n'
                          'Integration time: {}\n'
                          'Coadded spectra: {}\n'
                          'Wavelength [nm]\tIntensity [DN]'.format(self.spec_ctrl.int_time, self.spec_ctrl.coadd))

    def save_clear_raw(self, filename):
        """Save clear spectrum - not dark corrected"""
        # Sanity checks
        if self.doas_worker.wavelengths is None or self.doas_worker.clear_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.doas_worker.wavelengths) != len(self.doas_worker.clear_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        # Save file
        np.savetxt(filename, np.transpose([self.doas_worker.wavelengths, self.doas_worker.clear_spec_raw]),
                   header='Raw clear (Fraunhofer) spectrum\n'
                          '-Not dark-corrected\n'
                          'Integration time: {}\n'
                          'Coadded spectra: {}\n'
                          'Wavelength [nm]\tIntensity [DN]'.format(self.spec_ctrl.int_time, self.spec_ctrl.coadd))

    def save_plume_raw(self, filename):
        """Save plume spectrum - not dark corrected"""
        # Sanity checks
        if self.doas_worker.wavelengths is None or self.doas_worker.plume_spec_raw is None:
            raise ValueError('One or both attributes are NoneType. Cannot save.')
        if len(self.doas_worker.wavelengths) != len(self.doas_worker.plume_spec_raw):
            raise ValueError('Arrays are not the same length. Cannot save.')

        # Save file
        np.savetxt(filename, np.transpose([self.doas_worker.wavelengths, self.doas_worker.plume_spec_raw]),
                   header='Raw in-plume spectrum\n'
                          '-Not dark-corrected\n'
                          'Integration time: {}\n'
                          'Coadded spectra: {}\n'
                          'Wavelength [nm]\tIntensity [DN]'.format(self.spec_ctrl.int_time, self.spec_ctrl.coadd))

    def save_processed_spec(self, dark_path, clear_path, plume_path, doas_path):
        """Saves processed spectrum with all useful information"""
        # Sanity checks
        if self.doas_worker.wavelengths_cut is None or self.doas_worker.abs_spec_cut is None:
            return
        if self.doas_worker.ref_spec_types[0] not in self.doas_worker.ref_spec_fit:
            return

        # Save file
        np.savetxt(doas_path, np.transpose([self.doas_worker.wavelengths_cut, self.doas_worker.ref_spec_fit['SO2'],
                                                 self.doas_worker.abs_spec_cut]),
                   header='Processed DOAS spectrum\n'
                          'Dark spectrum: {}\nClear spectrum: {}\nPlume spectrum: {}\n'
                          'Shift: {}\nStretch: {}\n'
                          'Stray range [nm]: {}:{}\nFit window [nm]: {}:{}\n'
                          'Column density [ppm.m]: {}\n'
                          'Wavelength [nm]\tReference spectrum (fitted)\tAbsorbance spectrum'.format(
                          dark_path, clear_path, plume_path,
                          self.doas_worker.shift, self.doas_worker.stretch,
                          self.doas_worker.start_stray_wave, self.doas_worker.end_stray_wave,
                          self.doas_worker.start_fit_wave, self.doas_worker.end_fit_wave,
                          self.doas_worker.column_amount))

    def save_scan(self, save_dir):
        """Saves scan information"""
        try:
            # Save file
            np.savetxt(save_dir + self.scan_filename, np.transpose([self.scan_proc.scan_angles,
                                                                    self.scan_proc.column_densities]),
                       header='Processed scan details\n'
                              'Plume speed: {}\n'
                              'Plume distance: {}\n'
                              'Emission rate [kg/s]: {}\n'
                              'Emission rate [t/day]: {}\n'
                              'Scan angle [deg]\tColumn density [ppm.m]'.format(
                           self.scan_proc.plume_speed, self.scan_proc.plume_distance,
                           self.scan_proc.SO2_flux, self.scan_proc.flux_tons))
        except Exception as e:
            print(e)

