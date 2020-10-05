SpecScan

The SpecScan program is designed to perform DOAS retrievals of SO2 using an Ocean Optics spectrometer (which can be interfaced with via the seabreeze library) and a stepper motor.
Anaconda version 5.3 didn't work for me, but version 5.1.0 (5.2 may aso work) found here does: https://repo.anaconda.com/archive/

Requirements:
win32file: python -m pip install pywin32 or if using anaconda: conda install -c anaconda pywin32
seabreeze is best installed with Anaconda 'conda install -c poehlmann python-seabreeze'. 