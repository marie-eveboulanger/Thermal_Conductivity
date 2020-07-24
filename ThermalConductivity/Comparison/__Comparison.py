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
from ThermalConductivity.Utilities import Database as D
from ThermalConductivity import Visualization as V
from ThermalConductivity import Utilities as U

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
    __dict_measures = D.measurements_dict

    # Creation of a dictionnary to sort other info
    __dict_parameters = D.parameters_dict

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

        if filename is not None:
            filename = os.path.abspath(filename)
            data = U.read_file_treated(filename)
            header = U.read_header(filename)
            for key, value in data.items():
                self[key] = value
                self.measures.append(key)

            self["H"] = U.find_H(filename, header)
            self["date"] = U.find_date(filename, header)
            self["mount"] = U.find_mount(filename, header)
            self["sample"] = U.find_sample(filename, header)
            self["probe"] = U.find_probe(filename, header)

            self.parameters += ["H", "date", "mount", "sample", "probe"]

        for key, value in kwargs.items():
            self[key] = value
            if key not in self.parameters:
                self.parameters.append(key)

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
                if i != "Tp_Tm":
                    setattr(M, "__"+i, getattr(self, "__"+i)[key])
                else:
                    setattr(M, "__"+i, None)

            for i in self.parameters:
                setattr(M, "__"+i, getattr(self, "__"+i))

            M.measures = self.measures
            M.parameters = self.parameters
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
            self["dTx/T"] = self["dTx"]/self["T_av"]*100
        else:
            pass

        if "T_av" and "T0" and "dTx" in self.measures:
            self.measures.append("Resistance")
            self["Resistance"] = (self["T_av"]-self["T0"])/self["dTx"]
        else:
            pass

        if "dTx" and "dTy" in self.measures:
            self.measures.append("dTy/dTx")
            self["dTy/dTx"] = self["dTy"]/self["dTx"]*100
        else:
            pass

        if "kxx" and "kxy" in self.measures:
            self.measures.append("kxy/kxx")
            self["kxy/kxx"] = self["kxy"]/self["kxx"]*100
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

    # Creation of an internal dictionnary used to match measurements to their
    # respective axis titles to make the figures prettier.
    __dict_axis = V.axis_labels
    __dict_labels = V.legend_labels

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
        Used as a layer between the object and Visualization.Plot

        Parameters:
        ------------------------------------------------------------------------
        key:        string
                    The measurement to plot

        Kwargs:
        ------------------------------------------------------------------------
        show:       Bool
                    Determines if the figure is shown ore closed defaults to True

        parameters: list
                    list of parameters to be used for legends

        axis_fs:    Int
                    The axis labels fontsize

        fig:        matplotlib.figure
                    Used to draw on an existing figure, requires ax

        ax:         matplotlib.ax
                    Used to draw on an existing figure, requires fig

        x_axis:     string
                    The measurement to use as x-axis defaults to T_av
        """

        # Deal with kwargs
        # Tries to find sample name
        if "sample" in self.parameters:
            sample = self.measurements[0]["sample"]
            for m in self.measurements:
                if sample == m["sample"]:
                    pass
                else:
                    sample = None
                    break

        if "fig" in kwargs:
            return_fig = False
        else:
            return_fig = True

        if "x_axis" in kwargs:
            x_axis = kwargs["x_axis"]
            kwargs.pop("x_axis")
            if x_axis in self.measures:
                pass
            else:
                raise Exception("x_axis must be in self.measures")
        else:
            x_axis = "T_av"

        if "figtext" not in kwargs:
            kwargs["figtext"] = sample
        else:
            pass

        if "parameters" in kwargs:
            parameters = dict()
            parameters_list = kwargs["parameters"]
            kwargs.pop("parameters")
            parameters = []
            for p in parameters_list:
                if p in self.parameters:
                    for m in self.measurements:
                        params = dict()
                        params[p] = m[p]
                        parameters.append(params)
                else:
                    raise Exception("parameters must be in self.parameters")
        else:
            parameters = [dict() for i in self.measurements]

        kwargs["parameters"] = parameters

        if key != "Tp_Tm":

            if "show" in kwargs:
                show = kwargs["show"]
                kwargs["show"] = None
            else:
                kwargs["show"] = None
                show = True

            for i in range(len(self.measurements)):

                xdata, xkey = self.measurements[i][x_axis], x_axis
                ydata, ykey = self.measurements[i][key], key
                kwargs["parameters"] = parameters[i]

                if i == 0:
                    fig, ax = V.Plot(xdata, ydata, xkey, ykey, *args, **kwargs)
                elif i < len(self.measurements)-1:
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    V.Plot(xdata, ydata, xkey, ykey, *args, **kwargs)
                else:
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    kwargs["show"] = show
                    V.Plot(xdata, ydata, xkey, ykey, *args, **kwargs)
        else:
            if "show" in kwargs:
                show = kwargs["show"]
                kwargs["show"] = None
            else:
                kwargs["show"] = None
                show = True

            for i in range(len(self.measurements)):
                xdata, xkey = self.measurements[i][x_axis], x_axis
                ydata1, ykey1 = self.measurements[i]["Tp"], "Tp"
                ydata2, ykey2 = self.measurements[i]["Tm"], "Tm"
                kwargs["parameters"] = parameters[i]

                if i == 0:
                    kwargs["parameters"]["which"] = r"T$^{+}$"
                    fig, ax = V.Plot(xdata, ydata1, xkey,
                                     ykey1, *args, **kwargs)

                    kwargs["parameters"]["which"] = r"T$^{-}$"
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    fig, ax = V.Plot(xdata, ydata2, xkey,
                                     ykey2, *args, **kwargs)
                elif i < len(self.measurements)-1:
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    kwargs["parameters"]["which"] = r"T$^{+}$"
                    V.Plot(xdata, ydata1, xkey, ykey1, *args, **kwargs)

                    kwargs["parameters"]["which"] = r"T$^{-}$"
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    V.Plot(xdata, ydata2, xkey, ykey2, *args, **kwargs)

                else:
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    kwargs["parameters"]["which"] = r"T$^{+}$"
                    V.Plot(xdata, ydata1, xkey, ykey1, *args, **kwargs)

                    kwargs["show"] = show
                    kwargs["parameters"]["which"] = r"T$^{-}$"
                    kwargs["fig"], kwargs["ax"] = fig, ax
                    V.Plot(xdata, ydata2, xkey, ykey2, *args, **kwargs)

        if return_fig is False:
            return

        else:
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

        try:
            overwrite = kwargs["overwrite"]
            kwargs.pop("overwrite")
        except KeyError:
            overwrite = "ask"

        for key in measures:
            figures.append(self.Plot(key, *args, **kwargs)[0])

        if filename is not None:
            filename = os.path.abspath(filename)
            U.save_to_pdf(filename, figures, overwrite=overwrite)
        else:
            pass

        return

    def Plot_fancy(self, *args, **kwargs):
        """
        Just like Plot_all but with a fancy layout that is more suited to
        ipython notebooks
        """

        remove = ["T_av", "T0", "Tp", "Tm"]
        if len(self.measurements) > 1:
            remove.append("Tp_Tm")
        else:
            pass
        measures = [i for i in self.measures if i not in remove]
        ref_meas = ["kxx", "kxx/T", "kxy", "kxy/T", "dTx",
                    "dTx/T", "dTy/dTx", "Resistance"]
        measures = [i for i in ref_meas if i in measures]

        n = len(measures)
        fig, ax = V.create_grid(n)

        try:
            show = kwargs["show"]
            kwargs.pop("show")
        except KeyError:
            show = True

        # Tries to find sample name
        if "sample" in self.parameters:
            sample = self.measurements[0]["sample"]
            for m in self.measurements:
                if sample == m["sample"]:
                    pass
                else:
                    sample = None
                    break

        try:
            filename = kwargs["filename"]
            kwargs.pop("filename")
        except KeyError:
            filename = None

        try:
            overwrite = kwargs["overwrite"]
            kwargs.pop("overwrite")
        except KeyError:
            overwrite = "ask"

        for i in range(n):
            self.Plot(measures[i], *args, show=None,
                      fig=fig, ax=ax[i], **kwargs)

        if sample is not None:
            plt.suptitle(sample, y=0.95, fontsize=22)
        else:
            pass

        fig.tight_layout(rect=[0.01, 0.01, 1, 0.95])

        if filename is not None:
            filename = os.path.abspath(filename)
            U.save_to_pdf(filename, fig, overwrite=overwrite)
        else:
            pass

        if show is True:
            plt.show()
        elif show is False:
            plt.close()
        else:
            pass

        return fig, ax
