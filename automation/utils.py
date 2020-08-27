# -*- coding: utf-8 -*-

"""Utilities for pycam"""
import os
import numpy as np
from .setupclasses import SpecSpecs


def save_spectrum(wavelengths, spectrum, filename):
    """Saves spectrum as numpy .mat file
    wavelengths: NumPy array-like object
        Wavelength values held in array
    spectrum: NumPy array-like object
        Spectrum digital numbers held in array
    filename: str
        File path for saving
    """
    # Create lock file to secure file until saving is complete
    lock = filename.replace(SpecSpecs().file_ext, '.lock')
    open(lock, 'a').close()

    # Pack wavelengths and spectrum into single array
    spec_array = np.array([wavelengths, spectrum])

    # Save spectrum
    np.save(filename, spec_array)

    # Remove lock
    os.remove(lock)


def check_filename(filename, ext):
    """Checks filename to ensure it is as expected

    Parameters
    ----------
    filename: str
        full filename, expected to contain file extension <ext>
    ext: str
        expected filename extension to be checked
    """
    # Ensure filename is string
    if not isinstance(filename, str):
        raise ValueError('Filename must be in string format')

    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    # Split filename by .
    split_name = filename.split('.')

    # # Ensure filename contains exactly one . for file extension
    # if len(split_name) != 2:
    #     raise ValueError('Filename is not in the correct format. Name contained {} points'.format(len(split_name)-1))

    # Compare file extension to expected extension
    if split_name[-1] != ext:
        raise ValueError('Wrong file extension encountered')

    return


def write_file(filename, my_dict):
    """Writes all attributes of dictionary to file

    Parameters
    ----------
    filename: str
        file name to be written to
    my_dict: dict
        Dictionary of all data
    """
    # Check filename is legal
    try:
        check_filename(filename, 'txt')
    except ValueError:
        raise

    with open(filename, 'w') as f:
        # Loop through dictionary and write to file
        for key in my_dict:
            string = '{}={}\n'.format(key, my_dict[key])
            f.write(string)
    return


def read_file(filename, separator='=', ignore='#'):
    """Reads all lines of file separating into keys using the separator

        Parameters
        ----------
        filename: str
            file name to be written to
        separator: str
            string used to separate the key from its attribute
        ignore: str
            lines beginning with this string are ignored
            
        :returns
        data: dict
            dictionary of all attributes in file
    """
    # Check we are working with a text file
    check_filename(filename, 'txt')

    # Create empty dictionary to be filled
    data = dict()

    with open(filename, 'r') as f:

        # Loop through file line by line
        for line in f:

            # If line is start with ignore string then ignore line
            if line[0:len(ignore)] == ignore:
                continue

            try:
                # Split line into key and the key attribute
                key, attr = line.split(separator)[0:2]
            # ValueError will be thrown if nothing is after (or before) the equals sign. So we ignore these lines
            except ValueError:
                continue

            # Add attribute to dictionary, first removing any unwanted information at the end of the line
            # (including whitespace and #)
            data[key] = attr.split(ignore)[0].strip('\n').strip()

    return data


def format_time(time_obj, fmt):
    """Formats datetime object to string for use in filenames

    Parameters
    ----------
    time_obj: datetime.datetime
        Time to be converted to string"""
    return time_obj.strftime(fmt)
    # # Remove microseconds
    # time_obj = time_obj.replace(microsecond=0)
    #
    # # Return string format
    # return time_obj.isoformat().replace(':', '')
