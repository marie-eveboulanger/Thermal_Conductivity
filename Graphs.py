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


class Measurement():
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
    __dict_measures["kxx/T"] = ["kxx/T"]
    __dict_measures["Resistance"] = ["Resistance"]
    __dict_measures["dTx/T"] = ["dTx/T"]

    # Creation of a dictionnary to sort other info
    __dict_parameters = dict()
    __dict_parameters["H"] = ["H"]
    __dict_parameters["w"] = ["w"]
    __dict_parameters["t"] = ["t"]
    __dict_parameters["L"] = ["L"]
    __dict_parameters["Loc"] = ["BOT", "TOP", "Bot", "Top", "bot", "top"]
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
            setattr(self, "__H", H)
            self.parameters.append("H")

        if w is not None:
            setattr(self, "__w", w)
            self.parameters.append("w")

        if t is not None:
            setattr(self, "__t", t)
            self.parameters.append("t")

        if L is not None:
            setattr(self, "__L", L)
            self.parameters.append("L")

        if sample is not None:
            setattr(self, "__sample", sample)
            self.parameters.append("sample")

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
                            if len(self.measures) == 0:
                                raise Exception("No known measurements found")
                            else:
                                pass

            if hasattr(self, "__sample") is False:
                setattr(self, "__sample", "unknown")
            else:
                pass
            if hasattr(self, "__H") is False:
                setattr(self, "__H", "unknown")
            else:
                pass

        self.__add_measure()

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

    def __add_measure(self):
        if "T_av" and "kxx" in self.measures:
            self.measures.append("kxx/T")
            self["kxx/T"] = self["kxx"]/self["T_av"]
        else:
            pass

        if "T_av" and "dTx" in self.measures:
            self.measures.append("dTx/T")
            self["dTx/T"] = self["dTx"]/self["T_av"]
        else:
            pass

        if "T_av" and "T0" and "dTx" in self.measures:
            self.measures.append("Resistance")
            self["Resistance"] = (self["T_av"]-self["T0"])/self["dTx"]
        else:
            pass

        if "dTx" and "dTy" in self.measures:
            self.measures.append("dTy/dTx")
            self["dTy/dTx"] = self["dTy"]/self["dTx"]
        else:
            pass

        if "kxx" and "kxy" in self.measures:
            self.measures.append("kxy/kxx")
            self["kxy/kxx"] = self["kxy"]/self["kxx"]
        else:
            pass

        return


class Data_Set():
    """
    This class contains multiple Measurement objects and is used to compare them
    which means they must possess common attributes.
    """

    __list_measures = list()
    __list_parameters = list()
    __dict_axis = dict()
    __dict_axis["T_av"] = r"T ( K )"
    __dict_axis["T0"] = r"$T_0$ ( K )"
    __dict_axis["Tp"] = __dict_axis["T_av"]
    __dict_axis["Tm"] = __dict_axis["T_av"]
    __dict_axis["dTx"] = r"$\Delta T_{\rm x}$ ( K )"
    __dict_axis["kxx"] = r"$\kappa_{\rm xx}$ ( W / K m )"
    __dict_axis["kxx/T"] = r"$\kappa_{\rm xx}$/T ( W / K$^2$ m )"
    __dict_axis["dTx/T"] = r"$\Delta T_{\rm x} ( % )"
    __dict_axis["Resistance"] = r"(T-T$_0$)/$\Delta T_{\rm x}$"
    __dict_axis["kxy"] = r"$\kappa_{\rm xy}$ ( mW / K cm )"
    __dict_axis["kxy/kxx"] = r"$\kappa_{\rm xy}/\kappa_{\rm xx}$ ( % )"
    __dict_axis["dTy/dTx"] = r"$\Delta T_{\rm y}/\Delta T_{\rm x}$ ( % )"

    __dict_labels = dict()
    __dict_labels["H"] = r"H = %s T"
    __dict_labels["sample"] = r"Sample: %s"

    def __init__(self, measurements=None):
        """
        This class is meant to regroup multiple measurements to compare.

        Parameters:
        ------------------------------------------------------------------------

        measurements:   Measurement instance or list of instances
        """

        if measurements is None:
            measurements = []

        elif isinstance(measurements, Measurement):
            measurements = [measurements]

        elif isinstance(measurements, list):
            if len(measurements) == 0:
                pass
            else:
                for i in measurements:
                    if isinstance(i, Measurement) is False:
                        raise TypeError(
                            "All elements of the list must be "
                            "Measurements objects")
                    else:
                        pass
        else:
            raise TypeError(
                "Input must be None, a Measurement object "
                "or a list of Measurements")

        self.measurements = measurements
        if len(measurements) != 0:
            m = measurements[0]
            self.__list_measures = list(m._Measurement__dict_measures.keys())
            self.__list_parameters = list(
                m._Measurement__dict_parameters.keys())
        else:
            pass

        self.__find_measures()

        return

    def __add__(self, other):
        data = Data_Set()

        if isinstance(other, Data_Set) is True:
            data.Add_measurements(self.measurements)
            data.Add_measurements(other.measurements)
        else:
            raise TypeError("Can only add Data_Set objects")
        return data

    def __repr__(self):
        n = len(self.measurements)
        s = "Data set containing %i Measurements" % (n)
        return s

    def __getitem__(self, key):
        return self.measurements[key]

    def __setitem__(self, key, value):
        self.measurements[key] = value
        return

    def __delitem__(self, key):
        del self.measurements[key]
        return

    def Add_measurements(self, measurements):
        """
        Used to add measurements to the object after it is initialized

        Parameters:
        ------------------------------------------------------------------------

        measurements:   Measurement instance or list of instances
        """

        if measurements is None:
            measurements = []

        elif isinstance(measurements, Measurement):
            measurements = [measurements]

        elif isinstance(measurements, list):
            if len(measurements) == 0:
                pass
            else:
                for i in measurements:
                    if isinstance(i, Measurement) is False:
                        raise TypeError(
                            "All elements of the list must be "
                            "Measurements objects")
                    else:
                        pass
        else:
            raise TypeError(
                "Input must be None, a Measurement object "
                "or a list of Measurements")

        self.measurements += measurements
        if len(self.__list_measures) == 0:
            m = measurements[0]
            self.__list_measures = list(m._Measurement__dict_measures.keys())
            self.__list_parameters = list(
                m._Measurement__dict_parameters.keys())

        self.__find_measures()

        return

    def __find_measures(self):
        """
        Finds common measures in all Measurement objects, also does it for
        parameters
        """

        if hasattr(self, "measures") is False:
            measures = self.__list_measures
        else:
            measures = self.measures

        if hasattr(self, "parameters") is False:
            parameters = self.__list_parameters
        else:
            parameters = self.parameters

        if len(self.measurements) == 0:
            pass
        else:
            remove = []
            for m in self.measurements:
                for i in measures:
                    if i in m.measures:
                        pass
                    else:
                        remove.append(i)
            for i in remove:
                if i in measures:
                    measures.remove(i)
                else:
                    pass
            self.measures = measures

            remove = []
            for m in self.measurements:
                for i in parameters:
                    if i in m.parameters:
                        pass
                    else:
                        remove.append(i)
            for i in remove:
                if i in parameters:
                    parameters.remove(i)
                else:
                    pass
            self.parameters = parameters

        return

    def Plot(self, key, *args, **kwargs):
        """
        Plots data corresponding to key for all Measurement in Data_Set

        kwargs are all passed to ax.plot from matplotlib except the following:
        ------------------------------------------------------------------------
        show :  Bool
            Determines if the figure is shown ore closed defaults to True
        x_axis : str
            The key for the x-axis defaults to "T_av"
        parameter : str
            The key corresponding to the compared parameter defaults to "H"
        """

        # Looks for show as kwarg
        try:
            show = kwargs["show"]
            if type(show) is not bool:
                raise TypeError("show must be of type bool")
            else:
                kwargs.pop("show")
        except KeyError:
            show = True

        # Looks for x_axis as kwarg
        try:
            x_axis = kwargs["x_axis"]
            if x_axis not in self.measures:
                raise ValueError("x_axis must be in self.measures")
            else:
                kwargs.pop("x_axis")
        except KeyError:
            x_axis = "T_av"

        # Looks for parameter as kwarg
        try:
            parameter = kwargs["parameter"]
            if parameter not in self.parameters:
                raise ValueError("parameter must be in self.parameters")
            else:
                kwargs.pop("parameter")
        except KeyError:
            parameter = "H"

        if key in self.measures is False:
            raise ValueError("%s is not in self.measures") % (key)
        else:
            pass

        # Tries to find sample name
        if "sample" in self.parameters:
            sample = self.measurements[0]["sample"]
            for m in self.measurements:
                if sample == m["sample"]:
                    pass
                else:
                    sample = None
                    break

        fig, ax = plt.subplots()
        zero_line = 0

        # Draws the curves
        for m in self.measurements:
            ax.plot(m[x_axis], m[key], label=self.__dict_labels[parameter] % (
                m[parameter]), *args, **kwargs)
            if m[key].min()*m[key].max() < 0 and zero_line == 0:
                ax.plot(m[x_axis], 0*m[key], "--k", lw=2)
                zero_line += 1
            else:
                pass

        # If sample is the same for all measurements print it on the figure
        plt.figtext(0.1, 0.025, sample, fontsize=14)

        # Makes it pretty
        ax.set_xlabel(self.__dict_axis[x_axis], fontsize=16)
        ax.set_ylabel(self.__dict_axis[key], fontsize=16)
        ax.legend(fontsize=14)

        if show is False:
            plt.close()
        else:
            plt.show()
