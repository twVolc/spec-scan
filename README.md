SpecScan

The SpecScan program is designed to perform DOAS retrievals of SO2 using an Ocean Optics spectrometer (which can be interfaced with via the seabreeze library) and a stepper motor.
Anaconda version 5.3 didn't work for me, but version 5.1.0 (5.2 may aso work) found here does: https://repo.anaconda.com/archive/

Requirements:
python 3.8 (3.6 doesn't seem to work!)
win32file: python -m pip install pywin32 or if using anaconda: conda install -c anaconda pywin32
pyserial: python -m pip install pyserial
seabreeze is best installed with Anaconda 'conda install -c poehlmann python-seabreeze'
then: seabreeze_os_setup to install drivers


TODO - scaling for integration time, so this can be used for processing pycam spectrometer data - where integration time
TODO - changes through time