"""
This code is both a module and a python executable script meant to
analyse data from thermal conductivity experiment. This code should
support both the VTI and Tallahassee probes with and without magnetic
field. To use interactively import as module in Ipython or jupyter and
see docstrings provided within the code using function? in Ipython/juyter.
To use as a script simply run:
                python this_script.py /data/my_data.dat w t L
Where w, t and L are the dimensions of the sample. Using the script all
Figures will be saved as one pdf in /figures/my_data.pdf using the same
directory structure as the data folder.
"""
import sys
import os
import numpy as np
import numpy.polynomial.polynomial as npp
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from thermometry import seebeck_thermometry

################################################################################
#                          ____ _        _    ____ ____                        #
#                         / ___| |      / \  / ___/ ___|                       #
#                        | |   | |     / _ \ \___ \___ \                       #
#                        | |___| |___ / ___ \ ___) |__) |                      #
#                         \____|_____/_/   \_\____/____/                       #
################################################################################


class Conductivity():
    """
    This is the main class of the program. It contains all the
    data and other information about the sample. Also contains
    all the analysis functions for both probes.
    """

    # Creation of a dictionnary to sort data
    __dict_measures = dict()
    __dict_measures["T_av"] = ["T_av(K)", "Taverage(K)", "T (K)"]
    __dict_measures["Tp"] = ["T+(K)", "T+ (K)"]
    __dict_measures["Tm"] = ["T-(K)", "T- (K)"]
    __dict_measures["dTx"] = ["dTx(K)", "dTx (K)"]
    __dict_measures["kxx"] = ["kxx(W/Km)", "k_xx(W/Km)", "Kxx (W / K m)"]
    __dict_measures["kxy"] = ["kxy(W/mk)", "k_xy(W/Km)", "Kxy (W / K m)"]
    __dict_measures["dTy"] = ["dTy(K)", "dTy (K)"]
    __dict_measures["I"] = ["I(A)", "I (A)"]
    __dict_measures["dTabs"] = ["dTabs", "dT_abs"]
    __dict_measures["kxx/T"] = ["kxx/T"]
    __dict_measures["Resistance"] = ["Resistance"]
    __dict_measures["dTx/T"] = ["dTx/T"]
    __dict_measures["Tp_Tm"] = ["Tp_Tm"]

    # Creation of a dictionnary to sort raw data
    __dict_raw = dict()
    __dict_raw["T0"] = ["#T0(K)"]
    __dict_raw["I"] = ["I(A)"]
    __dict_raw["R+_0"] = ["R+_0(V)"]
    __dict_raw["R+_Q"] = ["R+_Q(V)"]
    __dict_raw["R-_0"] = ["R-_0(V)"]
    __dict_raw["R-_Q"] = ["R-_Q(V)"]
    __dict_raw["dTy_0"] = ["dTy_0(V)"]
    __dict_raw["dTy_Q"] = ["dTy_Q(V)"]
    __dict_raw["dTabs_0"] = ["Tabs_0(V)"]
    __dict_raw["dTabs_Q"] = ["Tabs_Q(V)"]
    __dict_raw["dTx_0"] = ["dTx_0(V)"]
    __dict_raw["dTx_Q"] = ["dTx_Q(V)"]

    # Creation of a dictionnary to sort other info
    __dict_parameters = dict()
    __dict_parameters["H"] = ["H"]
    __dict_parameters["w"] = ["w"]
    __dict_parameters["t"] = ["t"]
    __dict_parameters["L"] = ["L"]
    __dict_parameters["mount"] = ["BOT", "TOP", "Bot", "Top", "bot", "top"]
    __dict_parameters["sample"] = ["Sample", "sample"]
    __dict_parameters["date"] = ["Date", "date"]

    def __init__(self, filename=None, w=1e-6, t=1e-6, L=1e-6, sign=1):

        self.parameters = []
        dict_geo = {"w": w, "t": t, "L": L}
        for key, value in dict_geo.items():
            setattr(self, "__"+key, value)
            self.parameters.append(key)

        if sign in [1, -1]:
            self.sign = sign
        else:
            raise ValueError("Sign must be 1 or -1")

        self.__read_file(filename)

        if getattr(self, "__H") != "0.0":
            self.__symetrize()

        self.__analyze()

        return

    def __symetrize(self):
        sym = ["T0", "I", "R+_0", "R+_Q", "R-_0", "R-_Q",
               "dTabs_0", "dTabs_Q", "dTx_0", "dTx_Q"]
        anti_sym = ["dTy_0", "dTy_Q"]

        for i in sym:
            if i in self.raw_data:
                self[i] = 0.5*(self[i][0]+self[i][1])
            else:
                pass
        for i in anti_sym:
            if i in self.raw_data:
                self[i] = 0.5*(self[i][0]-self[i][1])
            else:
                pass
        return

    def __analyze(self):
        # Probe Tallahasse
        if self["probe"] == "Tallahasse":
            # Polyfit of R+ and R-
            Cp = npp.polyfit(np.log(self["R+_0"]), np.log(self["T0"]), 8)
            Cm = npp.polyfit(np.log(self["R-_0"]), np.log(self["T0"]), 8)
            index = np.where(self["R+_Q"] < self["R+_0"][-1])
            for i in self.raw_data:
                self[i] = np.delete(self[i], index)
            # Compute useful stuff
            self["Tp"] = np.exp(npp.polyval(np.log(self["R+_Q"]), Cp))
            self["Tm"] = np.exp(npp.polyval(np.log(self["R-_Q"]), Cm))
            self["dTx"] = self["Tp"]-self["Tm"]
            self["T_av"] = 0.5*(self["Tp"]+self["Tm"])
            alpha = self["w"]*self["t"]/self["L"]
            self["kxx"] = 5000*(self["I"])**2/self["dTx"]/alpha
            self.measures += ["T_av", "T0", "Tp", "Tm", "dTx", "kxx"]
            if self["H"] != "0.0":
                S = seebeck_thermometry((self["T_av"]+self["T0"])/2)
                self["dTy"] = self.sign*(self["dTy_Q"]-self["dTy_0"])/1000/S
                self["kxy"] = self["kxx"]*self["dTy"] / \
                    self["dTx"]*self["L"]/self["w"]
                self.measures += ["dTy", "kxy"]
        # VTI
        elif self["probe"] == "VTI":
            Q = 5000*self["I"]**2
            alpha = self["w"]*self["t"]/self["L"]
            T_av = 0*Q
            dT_abs = 0*Q
            dTx = 0*Q
            prev = T_av+1000
            # Loop to converge towards actual temperatures
            while abs(prev.sum()-T_av.sum()) > 1e-10:
                prev = T_av
                S1 = seebeck_thermometry((dT_abs/2+self["T0"]))
                S2 = seebeck_thermometry(self["T0"]+dT_abs+dTx/2)
                dT_abs = abs(self["dTabs_Q"]-self["dTabs_0"])/S1/1000
                dTx = (self["dTx_Q"]-self["dTx_0"])/S2/1000
                Tm = dT_abs+self["T0"]
                Tp = Tm+dTx
                T_av = Tm+dTx/2
            self["T_av"] = T_av
            self["Tp"] = Tp
            self["Tm"] = Tm
            self["dTx"] = dTx
            self["kxx"] = Q/self["dTx"]/alpha
            self.measures += ["T_av", "T0", "Tp", "Tm", "dTx", "kxx"]
            if self["H"] != "0.0":
                S = seebeck_thermometry(T_av)
                self["dTy"] = self.sign*(self["dTy_Q"]-self["dTy_0"])/S/1000
                self["kxy"] = self["kxx"]*self["dTy"] / \
                    self["dTx"]*self["L"]/self["w"]
                self.measures += ["dTy", "kxy"]

        return

    def __read_file(self, filename):
        """
        Used to read the file header and the data. Also detects the probe that
        has been used for the measurement. Also detects if the data is raw or
        already analyzed.
        """
        measures = []
        parameters = []
        raw_data = []
        # Converting filename to an absolute path if it is relative
        filename = os.path.abspath(filename)

        # Extracting info from filename
        # Date
        date = "-".join(filename.split("/")[-1].split(".")[-2].split("-")[-3:])
        setattr(self, "__date", date)
        parameters.append("date")

        # Magnetic field
        f = filename.split("/")[-1].split(".")
        H = ".".join([f[0][-2:].replace("-", ""), f[1][0]])
        setattr(self, "__H", H)
        parameters.append("H")

        # Sample name
        sample = list(filter(None, filename.split("/")[-1].split("-")))[-4]
        setattr(self, "__sample", sample)
        parameters.append("sample")

        # Mount
        mount = list(filter(None, filename.split("/")[-1].split("-")))[-5]
        setattr(self, "__mount", mount)
        parameters.append("mount")

        # Read the data
        if H == "0.0":
            data = np.genfromtxt(filename, delimiter="\t").T
        else:
            if len(filename.split("--")) == 1:
                filename2 = filename.replace("TS-", "TS--")
            else:
                filename2 = filename
                filename = filename2.replace("--", "-")

            data = np.genfromtxt(filename, delimiter="\t").T
            data2 = np.genfromtxt(filename2, delimiter="\t").T

        # Reading all the lines in the header
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

        # Makes sure there is a header
        if len(lines) == 0:
            raise Exception("No header detected, cannot analyze the data!")
        else:
            pass
        # Detect probe or already treated data
        for line in lines:
            l = list(filter(None, line.split("\t")))
            # Checks for comments shorter than the raw data header
            if len(l) < 6:
                for key, values in self.__dict_parameters.items():
                    if l[0] in values:
                        if hasattr(self, "__"+key) is False:
                            setattr(self, "__"+key, l[-1])
                            parameters.append(key)
                        else:
                            pass
                    else:
                        pass
            else:
                for key, values in self.__dict_raw.items():
                    for i in range(len(l)):
                        if l[i].strip() in values:
                            if H != "0.0":
                                setattr(self, "__"+key, [data[i], data2[i]])
                            else:
                                setattr(self, "__"+key, data[i])
                            raw_data.append(key)
                        else:
                            pass
                        if len(raw_data) == 0:
                            check_treated = True
                            self.__filetype = "treated"
                        else:
                            check_treated = False
                            self.__filetype = "raw_data"

                if check_treated is True:
                    for key, values in self.__dict_measures.items():
                        for i in range(len(l)):
                            if l[i].strip() in values:
                                setattr(self, "__"+key, data[i])
                                measures.append(key)
                            else:
                                pass
                            if len(measures) == 0:
                                raise Exception("No known measurements found")
                            else:
                                pass
                else:
                    pass

                if self.__filetype == "raw_data":
                    self.datatype = "Raw"
                    if "dTabs_0" in raw_data:
                        self["probe"] = "VTI"
                    else:
                        self["probe"] = "Tallahasse"
                elif self.__filetype == "treated":
                    self.datatype = "Treated"

                self.lines = lines
                self.parameters += parameters
                self.measures = measures
                self.raw_data = raw_data

        return

    def __getitem__(self, key):
        if type(key) is str:
            return getattr(self, "__"+key)
        else:
            M = Measurement()
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
