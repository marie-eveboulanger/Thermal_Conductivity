"""
This code is meant to be used interactively to compare data from different
samples. Ideally it would read the already analyzed data from a file that
contains the H, w,t and L values but they can also be specified manually
"""

################################################################################
#                 ___ __  __ ____   ___  ____ _____ ____                       #
#                |_ _|  \/  |  _ \ / _ \|  _ \_   _/ ___|                      #
#                 | || |\/| | |_) | | | | |_) || | \___ \                      #
#                 | || |  | |  __/| |_| |  _ < | |  ___) |                     #
#                |___|_|  |_|_|    \___/|_| \_\|_| |____/                      #
#                                                                              #
################################################################################

import sys
import os
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

################################################################################
#                 ____ _        _    ____ ____  _____ ____                     #
#                / ___| |      / \  / ___/ ___|| ____/ ___|                    #
#               | |   | |     / _ \ \___ \___ \|  _| \___ \                    #
#               | |___| |___ / ___ \ ___) |__) | |___ ___) |                   #
#                \____|_____/_/   \_\____/____/|_____|____/                    #
################################################################################


class Measurment():
    """
    This class is used to store the data from a single measurement in a standard
    format independantly of how to data file is structured.
    """

    # Creation of a dictionnary to sort data
    measures = dict()
    measures["T_av"] = ["T_av(K)", "Taverage(K)", "T (K)"]
    measures["T0"] = ["T0(K)", "T0 (K)"]
    measures["Tp"] = ["T+(K)", "T+ (K)"]
    measures["Tm"] = ["T-(K)", "T- (K)"]
    measures["dTx"] = ["dTx(K)", "dTx (K)"]
    measures["kxx"] = ["kxx(W/Km)", "k_xx(W/Km)", "Kxx (W / K m)"]
    measures["kxy"] = ["kxy(W/mk)", "k_xy(W/Km)", "Kxy (W / K m)"]
    measures["dTy"] = ["dTy(K)", "dTy (K)"]
    measures["I"] = ["I(A)", "I (A)"]
    measures["dTabs"] = ["dTabs", "dT_abs"]

    # Creation of a dictionnary to sort other info
    parameters = dict()
    parameters["H"] = ["H"]
    parameters["w"] = ["w"]
    parameters["t"] = ["t"]
    parameters["L"] = ["L"]
    parameters["sample"] = ["Sample", "sample"]

    def __init__(self, filename, H=None, w=None, t=None, L=None, sample=None):
        """
        Used to initialize a Measurement object

        Parameters:
        ------------------------------------------------------------------------
        filename :  str
                The path to the file containing the data to be read and stored.
                Can be relative or absolute.
        H, w, t, L :    int or float
                The value of the magnetic field used during the experiment and
                the geometric parameters of the sample.
        sample_name :   str
                The name of the sample.
        """

        # Sets the values using the info provided
        # if H is not None:
        #    self.H = H
        if w is not None:
            self.w = w
        if t is not None:
            self.t = t
        if L is not None:
            self.L = L
        if sample is not None:
            self.sample = sample

        filename = os.path.abspath(filename)
        raw_data = np.genfromtxt(filename, delimiter="\t", skip_header=1).T
        lines = []
        with open(filename) as f:
            for i, line in enumerate(f):
                l = line.split(" ")
                if l[0] != "#":
                    numbers = ["0", "1", "2", "3",
                               "4", "5", "6", "7", "8", "9"]
                    if l[0][0] in numbers:
                        break
                    else:
                        lines.append(l[0].strip())
                else:
                    lines.append(line.strip()[2:])

        #Should contain all the comment lines without trailing \n and starting #
        self.lines = lines

        # Separate misc comments from the actual header
        for line in self.lines:
            l = line.split("\t")
            if len(l) < 6:
                for key, values in self.parameters.items():
                    if l[0] in values:
                        if hasattr(self, key) is False:
                            setattr(self, key, l[-1])
            else:
                for key, values in self.measures.items():
                    for i in range(len(l)):
                        if l[i].strip() in values:
                            setattr(self, key, raw_data[i])

        if hasattr(self, "sample") is False:
            self.sample = "unknown"
        if hasattr(self, "H") is False:
            self.H = "unknown"

        return

    def __repr__(self):
        if self.H == "unknown":
            if self.sample != "unknown":
                string = "Measurement of %s at %s H" % (self.sample, self.H)
            else:
                string = "Measurement of %s sample at %s H" % (
                    self.sample, self.H)
        else:
            if self.sample != "unknown":
                string = "Measurement of %s at H=%sT" % (self.sample, self.H)
            else:
                string = "Measurement of %s sample at %s H" % (
                    self.sample, self.H)
        return string

    def __call__(self):
        if hasattr(self, "kxy") is False:
            data = np.array([self.T_av, self.T0, self.Tp,
                             self.Tm, self.dTx, self.kxx])
        else:
            data = np.array([self.T_av, self.T0, self.Tp,
                             self.Tm, self.dTx, self.kxx,self.dTy,self.kxy])

        return data

class Data_Set():
    """
    This class contains multiple Measurement objects and is used to compare them
    which means they must possess common attributes.
    """
