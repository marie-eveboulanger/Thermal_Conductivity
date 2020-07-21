"""
This module should contain technical stuff like readfile functions and date
handling that are required by other modules.
"""

import os
import datetime
import re
import numpy as np
import ThermalConductivity.Utilities.Database as D
from matplotlib.backends.backend_pdf import PdfPages


def get_symetric_file(filename, days=3):
    """
    Given a file name this function will try to find the file containg the data
    at the opposing magnetic field

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
    The name of the first file
    days:       int
    The number of days to check before or after in case measurement
    spanned multiple days
    """

    # Get the file and the directory
    filename = os.path.abspath(filename)
    directory, filename = tuple(os.path.split(filename))

    # Get all files in directory
    files = os.listdir(directory)

    # Remove all non-data files from list (logs, current functions etc)
    # by comparing the first letter in the filename since it starts with Data
    files = [f for f in files if f[0] == filename[0] and f != filename]

    # If only 0 file left
    if len(files) == 0:
        print("No symetric file found")
        filename2 = None

    # If 1 file left
    elif len(files) == 1:
        filename2 = files[0]

    # If there is 2 files checks for a treated file
    elif len(files) == 2:
        files = [f for f in files if f.find("treated") == -1]
        if len(files) == 1:
            filename2 = files[0]
        else:
            files.append("")

    # If there is more than 2 files use another method
    if len(files) > 2:
        H = find_H(os.path.join(directory, filename))
        if filename.find("--") == -1:
            filename2 = filename.replace(H, "-"+H)
        else:
            filename2 = filename.replace("-"+H, H)
        if_file = os.path.isfile(os.path.join(directory, filename2))
        if if_file is True:
            pass
        else:
            date = find_date(os.path.join(directory, filename))
            dates = generate_dates(date)
            for i in dates:
                filename3 = filename2.replace(date, i)
                if_file = os.path.isfile(
                    os.path.join(directory, filename3))
                if if_file is True:
                    filename2 = filename3
                    break
                else:
                    pass
            if if_file is True:
                pass
            else:
                print("No symetric file found")
                filename2 = None

    return os.path.join(directory, filename2)


def generate_dates(date, days=3):
    """
    Takes the date and generates dates for +/- n days. Usefull when looking for
    measurements that span multiple days.

    Parameters:
    ----------------------------------------------------------------------------
    date:   string
            Original date
    days:   int
            number of days to add/substract
    """
    # Convert to a useable format for datetime.date
    date = tuple([int(i) for i in date.split("-")])

    # Create datetime object
    date = datetime.date(*date)

    # Create range
    rnge = [i for i in range(-days, days+1) if i != 0]

    # Generate dates
    dates = [date+datetime.timedelta(i) for i in rnge]

    # Convert to string
    dates = [str(i) for i in dates]

    return dates


def find_H(filename, header=None):
    """
    Finds the magnetic field value in the filename or the file's header

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    header:     None or list of string
                Used to read a header that is already in memory, that way the
                file is only opened once which is more efficient
    """

    filename = os.path.abspath(filename)
    if header is None:
        header = read_header(filename)
    else:
        if type(header) is list:
            if type(header[0]) is str:
                pass
            else:
                raise TypeError(
                    "Header must be output of Utilities.read_header")
        else:
            raise TypeError("Header must be output of Utilities.read_header")
    if len(header) > 1:
        keys = D.parameters_dict["H"]
        H = None
        for h in header:
            if H is not None:
                break
            else:
                pass
            for k in keys:
                if h.find(k) != -1:
                    h = h.replace(k, "")
                    H = re.search(r"\d{1,2}\.\d{1}", h)
                    if H is not None:
                        H = H.group()
                        break
                    else:
                        pass
                else:
                    pass
    else:
        H = None

    if H is None:
        os.path.split(filename)[1]
        H = re.search(r"\d{1,2}\.\d{1}T", filename)
        if H is not None:
            H = H.group()[0:-1]
    else:
        pass

    return H


def find_date(filename, header=None):
    """
    Finds the date in the filename or the file's header

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    header:     None or list of string
                Used to read a header that is already in memory, that way the
                file is only opened once which is more efficient
    """

    filename = os.path.abspath(filename)
    if header is None:
        header = read_header(filename)
    else:
        if type(header) is list:
            if type(header[0]) is str:
                pass
            else:
                raise TypeError(
                    "Header must be output of Utilities.read_header")
        else:
            raise TypeError("Header must be output of Utilities.read_header")
    if len(header) > 1:
        keys = D.parameters_dict["date"]
        date = None
        for h in header:
            if date is not None:
                break
            else:
                pass
            for k in keys:
                if h.find(k) != -1:
                    h = h.replace(k, "")
                    date = re.search(r"\d{4}-\d{2}-\d{2}", h)
                    if date is not None:
                        date = date.group()
                        break
                    else:
                        pass
                else:
                    pass
    else:
        date = None

    if date is None:
        filename = os.path.split(filename)[1]
        date = re.search(r"\d{4}-\d{2}-\d{2}", filename)
        if date is not None:
            date = date.group()

    return date


def find_mount(filename, header=None):
    """
    Finds the mount in the filename or the file's header

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    header:     None or list of string
                Used to read a header that is already in memory, that way the
                file is only opened once which is more efficient
    """

    filename = os.path.abspath(filename)
    if header is None:
        header = read_header(filename)
    else:
        if type(header) is list:
            if type(header[0]) is str:
                pass
            else:
                raise TypeError(
                    "Header must be output of Utilities.read_header")
        else:
            raise TypeError("Header must be output of Utilities.read_header")
    if len(header) > 1:
        keys = D.parameters_dict["mount"]
        mount = None
        for h in header:
            if mount is not None:
                break
            else:
                pass
            for k in keys:
                if h.find(k) != -1:
                    h = h.replace(k, "")
                    mount = re.search(r"\w{3}", h)
                    if mount is not None:
                        mount = mount.group()
                        break
                    else:
                        pass
                else:
                    pass
    else:
        mount = None

    if mount is None:
        filename = os.path.split(filename)[1]
        mount = re.search(r"-\w{3}-", filename)
        if mount is not None:
            mount = mount.group()[1:-1]
    else:
        pass

    return mount


def find_sample(filename, header=None):
    """
    Finds the sample name in the filename or the file's header

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    header:     None or list of string
                Used to read a header that is already in memory, that way the
                file is only opened once which is more efficient
    """

    filename = os.path.abspath(filename)
    if header is None:
        header = read_header(filename)
    else:
        if type(header) is list:
            if type(header[0]) is str:
                pass
            else:
                raise TypeError(
                    "Header must be output of Utilities.read_header")
        else:
            raise TypeError("Header must be output of Utilities.read_header")
    if len(header) > 1:
        keys = D.parameters_dict["sample"]
        sample = None
        for h in header:
            if sample is not None:
                break
            else:
                pass
            for k in keys:
                if h.find(k) != -1:
                    h = h.replace(k, "")
                    sample = re.search(r"[\d\w]{4,}.*", h)
                    if sample is not None:
                        sample = sample.group().strip()
                        break
                    else:
                        pass
                else:
                    pass
    else:
        sample = None

    if sample is None:
        sample = os.path.split(filename)[1]
        sample = os.path.splitext(sample)[0]

        date = find_date(filename)
        H = find_H(filename)
        mount = find_mount(filename)

        sample = sample.replace(date, "")
        sample = sample.replace(H+"T", "")
        sample = sample.replace(mount, "")
        sample = sample.replace("Data", "")
        sample = re.search(r"[\d\w]{4,}", sample).group()
    else:
        pass

    return sample


def find_probe(filename, header):
    """
    Finds the probe used for the measurement

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    header:     None or list of string
                Used to read a header that is already in memory, that way the
                file is only opened once which is more efficient
    """

    if header is None:
        filename = os.path.abspath(filename)
        header = read_header(filename)
    else:
        pass

    header = list(filter(None, header[-1].strip().split("\t")))

    # Defaults to tallahassee changes only if VTI detected
    probe = "Tallahassee"
    for i in D.raw_data_dict["dTx_0"]:
        if i in header:
            probe = "VTI"
        else:
            pass

    return probe


def read_header(filename):
    """
    Reads the header of the specified file

    """
    filename = os.path.abspath(filename)

    with open(filename) as f:
        header = []
        for i, line in enumerate(f):
            if line[0] == "#":
                header.append(line[1:])
            elif line[0] != "#" and i == 0:
                header.append(line)
                break
            else:
                break

    return header


def read_file_raw(filename):
    """
    Reads the file and the file header and returns the data it found in a
    dictionary

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name and path of the file to be read
    """

    filename = os.path.abspath(filename)
    header = list(filter(None, read_header(filename)[-1].strip().split("\t")))
    data = np.genfromtxt(filename, delimiter="\t").T

    raw_data = dict()

    for key, values in D.raw_data_dict.items():
        for i in range(len(header)):
            if header[i] in values:
                raw_data[key] = data[i]
            else:
                pass

    return raw_data


def read_file_treated(filename):
    """
    Reads the file and the file header and returns the data it found in a
    dictionary

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name and path of the file to be read
    """

    filename = os.path.abspath(filename)
    header = list(filter(None, read_header(filename)[-1].strip().split("\t")))
    data = np.genfromtxt(filename, delimiter="\t").T

    measurements = dict()

    for key, values in D.raw_data_dict.items():
        for i in range(len(header)):
            if header[i] in values:
                measurements[key] = data[i]
            else:
                pass

    return measurements


def read_parameters(header):
    """
    Reads parameters stored in a header

    Parameters:
    ----------------------------------------------------------------------------
    header:     string
                Output from read_header
    """

    header = header[0:-1]
    parameters = dict()

    for key, values in D.parameters_dict.items():
        for h in header:
            line = h.strip().split("\t")
            param = line[0]
            value = line[-1]
            if param in values:
                parameters[key] = value
            else:
                pass

    return parameters


def write_to_file(filename, data, header, overwrite="ask", create_dir="ask"):
    """
    Writes data to a file and creates the needed directories.
    Also checks to prevent overwrites
    """

    # Gets the absolute path
    filename = os.path.abspath(filename)
    directory = os.path.split(filename)[0]

    # Checks if directory exists
    if os.path.isdir(directory) is True:
        pass
    elif os.path.isdir(directory) is False and create_dir == "ask":
        dir_ = directory
        answer = input(
            "Do you want to create the following directory: %s (Y/n)?" % dir_)
        if answer in ["Y", "y", "", "yes", "Yes"]:
            os.makedirs(directory)
            print("Created directory: %s" % directory)
        else:
            print("Could not create directory, aborting")
            return
    elif os.path.isdir(directory) is False and create_dir is True:
        os.makedirs(directory)
        print("Created directory: %s" % directory)

    elif os.path.isdir(directory) is False and create_dir is False:
        print("Could not create directory, aborting")
        return

    # Checks if file exists
    if os.path.isfile(filename) is False:
        np.savetxt(filename, data, delimiter="\t", header=header, fmt="%.6e")

    elif os.path.isfile(filename) is True and overwrite == "ask":
        answer = input(
            "Do you want to overwrite the following file: %s (Y/n)?" % filename)
        if answer in ["Y", "y", "", "yes", "Yes"]:
            np.savetxt(filename, data, delimiter="\t",
                       header=header, fmt="%.6e")
            print("Overwrote file: %s" % filename)
        else:
            print("Could not overwrite file, aborting")
            return
    elif os.path.isfile(filename) is True and overwrite is True:
        np.savetxt(filename, data, delimiter="\t", header=header, fmt="%.6e")

    elif os.path.isfile(filename) is True and overwrite is False:
        print("Could not overwrite file, aborting")
        return

    return


def save_to_pdf(filename, figures, overwrite="ask", create_dir="ask"):
    """
    Saves figure to a pdf and creates the needed directories.
    Also checks to prevent overwrites
    """

    if type(figures) is not list:
        figures = [figures]

    # Gets the absolute path
    filename = os.path.abspath(filename)
    directory = os.path.split(filename)[0]

    # Checks if directory exists
    if os.path.isdir(directory) is True:
        pass
    elif os.path.isdir(directory) is False and create_dir == "ask":
        dir_ = directory
        answer = input(
            "Do you want to create the following directory: %s (Y/n)?" % dir_)
        if answer in ["Y", "y", "", "yes", "Yes"]:
            os.makedirs(directory)
            print("Created directory: %s" % directory)
        else:
            print("Could not create directory, aborting")
            return
    elif os.path.isdir(directory) is False and create_dir is True:
        os.makedirs(directory)
        print("Created directory: %s" % directory)

    elif os.path.isdir(directory) is False and create_dir is False:
        print("Could not create directory, aborting")
        return

    # Checks if file exists
    if os.path.isfile(filename) is False:
        pp = PdfPages(filename)
        for fig in figures:
            pp.savefig(fig)
        pp.close()

    elif os.path.isfile(filename) is True and overwrite == "ask":
        answer = input(
            "Do you want to overwrite the following file: %s (Y/n)?" % filename)
        if answer in ["Y", "y", "", "yes", "Yes"]:
            pp = PdfPages(filename)
            for fig in figures:
                pp.savefig(fig)
            pp.close()
            print("Overwrote file: %s" % filename)
        else:
            print("Could not overwrite file, aborting")
            return
    elif os.path.isfile(filename) is True and overwrite is True:
        pp = PdfPages(filename)
        for fig in figures:
            pp.savefig(fig)
        pp.close()

    elif os.path.isfile(filename) is True and overwrite is False:
        print("Could not overwrite file, aborting")
        return

    return
