import sys
import os
import datetime
import numpy as np
import numpy.polynomial.polynomial as npp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ThermalConductivity.Analysis import Functions as F
from ThermalConductivity import Utilities as U
from ThermalConductivity.Utilities import Database as D
from ThermalConductivity.Thermometry import seebeck_thermometry
from ThermalConductivity import Visualization as V

################################################################################
#                          ____ _        _    ____ ____                        #
#                         / ___| |      / \  / ___/ ___|                       #
#                        | |   | |     / _ \ \___ \___ \                       #
#                        | |___| |___ / ___ \ ___) |__) |                      #
#                         \____|_____/_/   \_\____/____/                       #
################################################################################


class Log():
    """
    This class is meant to read data from a log file for debugging
    """

    def __init__(self, filename):

        # Importing dictionaries
        self["dict_measures"] = D.log_data_dict
        self["dict_parameters"] = D.parameters_dict

        # Initializing lists
        self.measures = []
        self.parameters = []

        # Read the file
        filename = os.path.abspath(filename)
        header = U.read_header(filename)

        # Find the parameters
        self["H"] = U.find_H(filename)
        self["date"] = U.find_date(filename)
        self["mount"] = U.find_mount(filename)
        self["probe"] = U.find_probe(filename)
        self["sample"] = U.find_sample(filename, header)
        self.parameters.append(["H", "date", "mount", "probe", "sample"])

        # Find the measurements
        log_data = U.read_file_log(filename)

        for key, value in log_data:
            self[key] = value
            self.measures.append(key)

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
            x_axis = "Time"

        if "figtext" not in kwargs:
            kwargs["figtext"] = self["sample"]
        else:
            pass

        if "parameters" in kwargs:
            parameters = dict()
            parameters_list = kwargs["parameters"]
            kwargs.pop("parameters")
            for p in parameters_list:
                if p in self.parameters:
                    parameters[p] = self[p]
                else:
                    raise Exception("parameters must be in self.parameters")
        else:
            parameters = dict()

        kwargs["parameters"] = parameters

        # Plot the data
        xdata, xkey = self[x_axis], x_axis
        ydata, ykey = self[key], key

        fig, ax = V.Plot(xdata, ydata, xkey, ykey, *args, **kwargs)

        if return_fig is False:
            return
        else:
            return fig, ax
