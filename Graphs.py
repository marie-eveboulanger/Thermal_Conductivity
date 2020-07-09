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
    __dict_measures = dict()
    __dict_measures["T_av"] = ["T_av(K)", "Taverage(K)", "T (K)"]
    __dict_measures["T0"] = ["T0(K)", "T0 (K)"]
    __dict_measures["Tp"] = ["T+(K)", "T+ (K)"]
    __dict_measures["Tm"] = ["T-(K)", "T- (K)"]
    __dict_measures["dTx"] = ["dTx(K)", "dTx (K)"]
    __dict_measures["kxx"] = ["kxx(W/Km)", "k_xx(W/Km)", "Kxx (W / K m)"]
    __dict_measures["kxy"] = ["kxy(W/mk)", "k_xy(W/Km)", "Kxy (W / K m)"]
    __dict_measures["dTy"] = ["dTy(K)", "dTy (K)"]
    __dict_measures["I"] = ["I(A)", "I (A)"]
    __dict_measures["dTabs"] = ["dTabs", "dT_abs"]

    # Creation of a dictionnary to sort other info
    __dict_parameters = dict()
    __dict_parameters["H"] = ["H"]
    __dict_parameters["w"] = ["w"]
    __dict_parameters["t"] = ["t"]
    __dict_parameters["L"] = ["L"]
    __dict_parameters["sample"] = ["Sample", "sample"]

    def __init__(self, filename=None, H=None, w=None, t=None,
                 L=None, sample=None):
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

        self.measures = []
        self.parameters = []
        # Sets the values using the info provided
        if H is not None:
            self.H = H

        if w is not None:
            self.w = w

        if t is not None:
            self.t = t

        if L is not None:
            self.L = L

        if sample is not None:
            self.sample = sample

        if filename is not None:
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

            # Should contain all the comment lines without trailing \n and
            #starting #
            self.lines = lines

            # Separate misc comments from the actual header
            for line in self.lines:
                l = line.split("\t")
                if len(l) < 6:
                    for key, values in self.__dict_parameters.items():
                        if l[0] in values:
                            if hasattr(self, key) is False:
                                setattr(self, "__"+key, l[-1])
                                self.parameters.append(key)
                            else:
                                pass
                        else:
                            pass
                else:
                    for key, values in self.__dict_measures.items():
                        for i in range(len(l)):
                            if l[i].strip() in values:
                                setattr(self, "__"+key, raw_data[i])
                                self.measures.append(key)
                            else:
                                pass

            if hasattr(self, "__sample") is False:
                self.__sample = "unknown"
            else:
                pass
            if hasattr(self, "__H") is False:
                self.__H = "unknown"
            else:
                pass

        return

    def __repr__(self):
        H = getattr(self, "__H")
        S = getattr(self, "__sample")
        if H == "unknown":
            if S != "unknown":
                string = "Measurement of %s at %s H" % (S, H)
            else:
                string = "Measurement of %s sample at %s H" % (S, H)
        else:
            if S != "unknown":
                string = "Measurement of %s at H=%sT" % (S, H)
            else:
                string = "Measurement of %s sample at %s H" % (S, H)

        return string

    def __call__(self):
        """
        Returns the transposed raw data
        """
        data = np.array([getattr(self, "__"+key) for key in self.measures])
        return data

    def __getitem__(self, key):
        if type(key) is str:
            return getattr(self, "__"+key)
        else:
            M = Measurment()
            for i in self.measures:
                setattr(M, "__"+i, getattr(self, "__"+i)[key])

            for i in self.parameters:
                setattr(M, "__"+i, getattr(self, "__"+i))

            M["measures"] = self.measures
            M["parameters"] = self.parameters
            return M

    def __setitem__(self, key, value):
        if type(key) is str:
            setattr(self, "__"+key, value)
        else:
            pass
        return

    def __delitem__(self, key):
        delattr(self, "__"+key)
        return


class Data_Set():
    """
    This class contains multiple Measurement objects and is used to compare them
    which means they must possess common attributes.
    """
