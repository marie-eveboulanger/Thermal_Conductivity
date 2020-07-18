"""
This module should contain technical stuff like readfile functions and date
handling that are required by other modules.
"""

import os
import datetime
import re


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
        H = find_H(filename)
        if filename.find("--") == -1:
            filename2 = filename.replace(H, "-"+H)
        else:
            filename2 = filename.replace("-"+H, H)
        if_file = os.path.isfile(os.path.join(directory, filename2))
        if if_file is True:
            pass
        else:
            date = find_date(filename2)
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

    return os.path.join(directory,filename2)


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


def find_H(filename):
    """
    Finds the magnetic field value in the filename

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    """

    os.path.split(filename)[1]
    H = re.search(r"\d{1,2}\.\d{1}T", filename)
    if H is not None:
        H = H.group()[0:-1]

    return H


def find_date(filename):
    """
    Finds the date in the filename

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    """

    filename = os.path.split(filename)[1]
    date = re.search(r"\d{4}-\d{2}-\d{2}", filename)
    if date is not None:
        date = date.group()

    return date


def find_mount(filename):
    """
    Finds the mount in the filename

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    """

    filename = os.path.split(filename)[1]
    mount = re.search(r"-\w{3}-", filename)
    if mount is not None:
        mount = mount.group()[1:-1]

    return mount


def find_sample(filename):
    """
    Finds the sample name in the filename

    Parameters:
    ----------------------------------------------------------------------------
    filename:   string
                The name of the file to search
    """

    filename = os.path.split(filename)[1]
    sample = os.path.splitext(filename)[0]

    date = find_date(filename)
    H = find_H(filename)
    mount = find_mount(filename)

    sample = sample.replace(date, "")
    sample = sample.replace(H+"T", "")
    sample = sample.replace(mount, "")
    sample = sample.replace("Data", "")
    sample = re.search(r"[\d\w]{4,}", sample).group()

    return sample
