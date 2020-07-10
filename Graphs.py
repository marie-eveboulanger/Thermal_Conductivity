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
    __dict_measures["Tp_Tm"] = ["Tp_Tm"]

    # Creation of a dictionnary to sort other info
    __dict_parameters = dict()
    __dict_parameters["H"] = ["H"]
    __dict_parameters["w"] = ["w"]
    __dict_parameters["t"] = ["t"]
    __dict_parameters["L"] = ["L"]
    __dict_parameters["Loc"] = ["BOT", "TOP", "Bot", "Top", "bot", "top"]
    __dict_parameters["sample"] = ["Sample", "sample"]
    __dict_parameters["date"] = ["Date", "date"]

    def __init__(self, filename=None, **kwargs):
        """
        Used to initialize a Measurement object

        Parameters:
        ------------------------------------------------------------------------
        filename :  str
                The path to the file containing the data to be read and stored.
                Can be relative or absolute.

        Useful kwargs:
        ------------------------------------------------------------------------
        H, w, t, L :    int or float
                The value of the magnetic field used during the experiment and
                the geometric parameters of the sample. Note that kwargs are
                case sensitive and are used to populate the parameters attribute
                of the object
        sample :   str
                The name of the sample.
        """

        self.measures = []
        self.parameters = []
        # Sets the values using the info provided
        for key, value in kwargs.items():
            setattr(self, "__"+key, value)
            self.parameters.append(key)
            try:
                self.__dict_parameters[key]
            except KeyError:
                self.__dict_parameters[key] = [key]

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
                            if hasattr(self, "__"+key) is False:
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

        if "Tp" and "Tm" in self.measures:
            self.measures.append("Tp_Tm")
            self["Tp_Tm"] = None

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
    __dict_axis["dTx/T"] = r"$\Delta T_{\rm x}$/T ( % )"
    __dict_axis["Resistance"] = r"(T-T$_0$)/$\Delta T_{\rm x}$"
    __dict_axis["kxy"] = r"$\kappa_{\rm xy}$ ( mW / K cm )"
    __dict_axis["kxy/kxx"] = r"$\kappa_{\rm xy}/\kappa_{\rm xx}$ ( % )"
    __dict_axis["dTy/dTx"] = r"$\Delta T_{\rm y}/\Delta T_{\rm x}$ ( % )"
    __dict_axis["Tp_Tm"] = __dict_axis["T_av"]

    __dict_labels = dict()
    __dict_labels["H"] = r"H = %sT"
    __dict_labels["sample"] = r"Sample: %s"
    __dict_labels["date"] = r"%s"

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
                if show is not None:
                    raise TypeError("show must be of type bool or None")
                else:
                    pass
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
            parameters = kwargs["parameters"]
            if type(parameters) is not list:
                if type(parameters) is str:
                    parameters = [parameters]
                else:
                    raise TypeError("Parameter must be a string")
            else:
                pass
            for parameter in parameters:
                if parameter not in self.parameters:
                    raise ValueError("parameter must be in self.parameters")
                else:
                    pass
            kwargs.pop("parameters")

        except KeyError:
            parameters = ["H"]

        # Sets the label according to parameters
        labels = []
        for parameter in parameters:
            try:
                label = self.__dict_labels[parameter]
            except KeyError:
                label = "%s"
            labels.append(label)
        label = ", ".join(labels)

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
        if key == "Tp_Tm":
            for m in self.measurements:
                p = [m[parameter] for parameter in parameters]
                ax.plot(m[x_axis], m["Tp"], label="T+ "+label %
                        tuple(p), *args, **kwargs)
                ax.plot(m[x_axis], m["Tm"], label="T- "+label %
                        tuple(p), *args, **kwargs)
        else:
            for m in self.measurements:
                p = [m[parameter] for parameter in parameters]
                x_data = m[x_axis]
                if key in ["kxy", "kxy/T"]:
                    y_data = 10*m[key]
                else:
                    y_data = m[key]
                ax.plot(m[x_axis], m[key], label=label %
                        tuple(p), *args, **kwargs)
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
        elif show is True:
            plt.show()

        return fig, ax

    def Plot_all(self, *args, **kwargs):
        """
        Plots all non trivial measures, all the same kwargs as Data_Set.Plot
        with the addition of filename to save the file.
        """

        remove = ["T_av", "T0", "Tp", "Tm"]
        if len(self.measurements) > 1:
            remove.append("Tp_Tm")
        else:
            pass
        measures = [i for i in self.measures if i not in remove]
        figures = []

        try:
            filename = kwargs["filename"]
            kwargs.pop("filename")
        except KeyError:
            filename = None

        for key in measures:
            figures.append(self.Plot(key, *args, **kwargs)[0])

        if filename is not None:
            filename = os.path.abspath(filename)
            pp = PdfPages(filename)
            for i in figures:
                pp.savefig(i)
            pp.close()
        else:
            pass

        return
