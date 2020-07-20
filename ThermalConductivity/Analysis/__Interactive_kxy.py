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
import datetime
import numpy as np
import numpy.polynomial.polynomial as npp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ThermalConductivity.Analysis import Functions as F
from ThermalConductivity import Utilities as U
from ThermalConductivity.Utilities import Database as D
from ThermalConductivity.Thermometry import seebeck_thermometry

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
    __dict_measures = D.measurements_dict

    # Creation of a dictionnary to sort raw data
    __dict_raw = D.raw_data_dict

    # Creation of a dictionnary to sort other info
    __dict_parameters = D.parameters_dict

    # Creation of an internal dictionnary used to match measurements to their
    # respective axis titles to make the figures prettier.
    __list_measures = list()
    __list_parameters = list()
    __dict_axis = dict()
    __dict_axis["T_av"] = r"T ( K )"
    __dict_axis["T0"] = r"$T_0$ ( K )"
    __dict_axis["Tp"] = __dict_axis["T_av"]
    __dict_axis["Tm"] = __dict_axis["T_av"]
    __dict_axis["kxx"] = r"$\kappa_{\rm xx}$ ( W / K m )"
    __dict_axis["dTx"] = r"$\Delta T_{\rm x}$ ( K )"
    __dict_axis["dTy"] = r"$\Delta T_{\rm y}$ ( K )"
    __dict_axis["kxx/T"] = r"$\kappa_{\rm xx}$/T ( W / K$^2$ m )"
    __dict_axis["dTx/T"] = r"$\Delta T_{\rm x}$/T ( % )"
    __dict_axis["Resistance"] = r"(T-T$_0$)/$\Delta T_{\rm x}$"
    __dict_axis["kxy"] = r"$\kappa_{\rm xy}$ ( mW / K cm )"
    __dict_axis["kxy/kxx"] = r"$\kappa_{\rm xy}/\kappa_{\rm xx}$ ( % )"
    __dict_axis["dTy/dTx"] = r"$\Delta T_{\rm y}/\Delta T_{\rm x}$ ( % )"
    __dict_axis["Tp_Tm"] = __dict_axis["T_av"]
    __dict_axis["T0_fit"] = __dict_axis["T0"]
    __dict_axis["I_fit"] = r"I ( mA )"

    # Same principle then before but for curve labels
    __dict_labels = dict()
    __dict_labels["H"] = r"H = %sT"
    __dict_labels["sample"] = r"Sample: %s"
    __dict_labels["date"] = r"%s"

    def __init__(self, filename=None, w=1e-6, t=1e-6, L=1e-6, sign=1, **kwargs):

        # Check for some specific kwargs
        try:
            self["force_kxy"] = kwargs["force_kxy"]
            kwargs.pop("force_kxy")
            if type(self["force_kxy"]) is bool:
                pass
            else:
                raise TypeError("force_kxy must be True or False")
        except KeyError:
            self["force_kxy"] = False

        try:
            self["symmetrize"] = kwargs["symmetrize"]
            kwargs.pop("symmetrize")
            if type(self["symmetrize"]) is bool:
                pass
            else:
                raise TypeError("symmetrize must be True or False")
        except KeyError:
            self["symmetrize"] = True

        if sign in [1, -1]:
            self["sign"] = sign
        else:
            raise ValueError("Sign must be 1 or -1")

        if filename is not None:
            filename = os.path.abspath(filename)
            header = U.read_header(filename)

            # Find info
            self["H"] = U.find_H(filename, header)
            self["date"] = U.find_date(filename, header)
            self["mount"] = U.find_mount(filename, header)
            self["sample"] = U.find_sample(filename, header)
            self["probe"] = U.find_probe(filename, header)
            self["w"] = w
            self["t"] = t
            self["L"] = L

            self.parameters = ["H", "date", "mount",
                               "sample", "probe", "w", "t", "L"]

            # If symetrize is True
            if self["H"] != "0.0" and self["symmetrize"] is True:
                filename2 = U.get_symetric_file(filename)
                raw_data = self.__Symmetrize(filename, filename2)

            else:
                raw_data = U.read_file_raw(filename)

            for key, values in raw_data.items():
                self[key] = values


            self.measures = []
            self.__Analyze()
            self.__add_measure()

        # Remaining kwargs are set as parameters
        for key, value in kwargs.items():
            self[key] = value
            self.parameters.append(key)

        return

    def __Symmetrize(self, filename, filename2):
        anti_sym = ["dTy_0", "dTy_Q"]

        if filename.find("--") != -1:
            filename, filename2 = filename2, filename
        else:
            pass

        data = U.read_file_raw(filename)
        data2 = U.read_file_raw(filename2)

        sym_data = dict()

        for key, values in data:
            if key in anti_sym:
                sym_data[key] = 0.5*(data[key]-data2[key])
            else:
                sym_data[key] = 0.5*(data[key]+data[key])

        return sym_data

    def __Analyze(self):
        # Probe Tallahasse
        if self["probe"] == "Tallahasse":
            # Cut the uncalibrated points
            index = np.where(self["R+_Q"] < self["R+_0"][-1])
            for i in self.raw_data:
                self[i] = np.delete(self[i], index)

            # Compute useful stuff
            # Get I and T0
            I = self["I"]
            T0 = self["T0"]

            # Compute T+ and T-
            Tp = F.tallahassee_temp(self["R+_0"], self["R+_Q"], T0)
            Tm = F.tallahassee_temp(self["R-_0"], self["R-_Q"], T0)

            # Compute T_av dTx and kxx
            T_av = 0.5*(Tp+Tm)
            dTx = (Tp-Tm)
            kxx = F.compute_kxx(I, dTx, self["w"], self["t"], self["L"])

            # Store values in self
            self["kxx"] = kxx
            self["dTx"] = dTx
            self["T_av"] = T_av
            self["Tp"] = Tp
            self["Tm"] = Tm
            self.measures += ["T_av", "T0", "Tp", "Tm", "dTx", "kxx"]

            # Compute the transverse stuff
            if self["H"] != "0.0" or self["force_kxy"] is True:
                # Compute dTy
                Tr = T0+T_av/2  # Reference tempereature for the thermocouple
                dTy = F.compute_thermocouple(self["dTy_0"], self["dTy_Q"], Tr)
                dTy *= self["sign"]  # Apply the sign

                # Compute kxy
                kxy = F.compute_kxy(kxx, dTx, dTy, self["w"], self["L"])

                # Store in self
                self["dTy"] = dTy
                self["kxy"] = kxy
                self.measures += ["dTy", "kxy"]

        # VTI
        elif self["probe"] == "VTI":

            # Importing data
            dTabs_0, dTabs_Q = self["dTabs_0"], self["dTabs_Q"]
            dTx_0, dTx_Q = self["dTx_0"], self["dTx_Q"]
            T0 = self["T0"]
            I = self["I"]

            # Computing everything
            result = F.vti_calibration_loop(dTabs_0, dTabs_Q, dTx_0, dTx_Q, T0)
            kxx = F.compute_kxx(I, result[1], self["w"], self["t"], self["L"])

            # Storing in self
            self["kxx"] = kxx
            self["T_av"] = result[0]
            self["dTx"] = result[1]
            self["Tp"] = result[2]
            self["Tm"] = result[3]
            self.measures += ["T_av", "T0", "Tp", "Tm", "dTx", "kxx"]

            if self["H"] != "0.0" or self["force_kxy"] is True:
                # Compute dTy
                Tr = (T0+self["T_av"])/2  # Reference temp for the thermocouple
                dTy = F.compute_thermocouple(self["dTy_0"], self["dTy_Q"], Tr)
                dTy *= self["sign"]  # Apply the sign

                # Compute kxy
                kxy = F.compute_kxy(
                    kxx, self["dTx"], dTy, self["w"], self["L"])

                # Store in self
                self["dTy"] = dTy
                self["kxy"] = kxy
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
        self["filename"] = filename

        # Extracting info from filename
        # Date
        date = "-".join(filename.split("/")[-1].split(".")[-2].split("-")[-3:])
        setattr(self, "__date", date)
        parameters.append("date")

        # Magnetic field
        f = filename.split("/")[-1].split(".")
        H = ".".join([f[0][-2:].replace("-", ""), f[1][0]])
        if "H" != "0.0" and self["symmetrize"] is False:
            if len(filename.split("--")) == 1:
                pass
            else:
                H = "-"+H
        setattr(self, "__H", H)
        parameters.append("H")
        if H == "0.0":
            self["symmetrize"] = False
        else:
            pass

        # Sample name
        if hasattr(self, "__sample") is False:
            sample = list(filter(None, filename.split("/")[-1].split("-")))[-4]
            setattr(self, "__sample", sample)
            parameters.append("sample")
        else:
            pass

        # Mount
        mount = list(filter(None, filename.split("/")[-1].split("-")))[-5]
        setattr(self, "__mount", mount)
        parameters.append("mount")

        # Read the data
        if H == "0.0" or self["symmetrize"] is False:
            data = np.genfromtxt(filename, delimiter="\t").T
        else:
            if len(filename.split("--")) == 1:
                filename2 = filename.replace("TS-", "TS--")
                exist = os.path.isfile(filename2)
                if exist is True:
                    pass
                else:
                    filename3 = filename2
                    dates = self.__Dates(date)
                    for i in dates:
                        filename3 = filename2.replace(date, i)
                        exist = os.path.isfile(filename3)
                        if exist is True:
                            filename2 = filename3
                            break
                        else:
                            pass
                    if exist is False:
                        filename2 = filename
                    else:
                        filename2 = filename3
            else:
                filename2 = filename
                filename = filename2.replace("--", "-")
                exist = os.path.isfile(filename)
                if exist is True:
                    pass
                else:
                    filename3 = filename
                    dates = self.__Dates(date)
                    for i in dates:
                        filename3 = filename.replace(date, i)
                        exist = os.path.isfile(filename3)
                        if exist is True:
                            filename = filename3
                            break
                        else:
                            pass
                    if exist is False:
                        filename = filename2
                    else:
                        filename = filename3

            if filename == filename2:
                self["symmetrize"] = False
            else:
                pass
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
                            if H != "0.0" and self["symmetrize"] is True:
                                setattr(self, "__"+key, [data[i], data2[i]])
                            else:
                                setattr(self, "__"+key, data[i])
                            raw_data.append(key)
                        else:
                            pass
                        if len(raw_data) == 0:
                            check_treated = True
                            self["filetype"] = "treated"
                        else:
                            check_treated = False
                            self["filetype"] = "raw_data"

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

                if self["filetype"] == "raw_data":
                    self.datatype = "Raw"
                    if "dTabs_0" in raw_data:
                        self["probe"] = "VTI"
                    else:
                        self["probe"] = "Tallahasse"
                elif self["filetype"] == "treated":
                    self.datatype = "Treated"

                self.lines = lines
                self.parameters += parameters
                self.measures = measures
                self.raw_data = raw_data

        return

    def __Dates(self, date):
        d = datetime.date(*tuple([int(i) for i in date.split("-")]))
        dates = ["%s" % (d+datetime.timedelta(i))
                 for i in range(-2, 3) if i != 0]
        return dates

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

    def __create_grid(self, measures):
        n = len(measures)

        if n % 2 == 0:
            N = n//2
            fig, ax = plt.subplots(N, 2, figsize=(16, N*4.5))
            axes = ax.flatten().tolist()
        else:
            N = n//2+1
            s = (N, 4)
            fig, ax = plt.subplots(N, 2, figsize=(16, N*4.5))
            axes = []
            loc = (0, 0)
            for i in range(n):

                if i != 0:
                    if i != n-1:
                        if int(loc[1]) == 0:
                            loc = (loc[0], 2)
                        elif int(loc[1]) == 2:
                            loc = (loc[0]+1, 0)
                    else:
                        loc = (N-1, 1)
                else:
                    pass
                axes.append(plt.subplot2grid(s, loc, colspan=2, fig=fig))

        return fig, axes

    def Plot(self, key, *args, **kwargs):
        """
        Plots data corresponding to key.

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
                    kwargs.pop("show")
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

        # Looks for axis_fontsize as kwarg
        try:
            axis_fs = kwargs["axis_fontsize"]
            kwargs.pop("axis_fontsize")
        except KeyError:
            axis_fs = 16

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
            parameters = []

        # Sets the label according to parameters
        labels = []
        for parameter in parameters:
            try:
                label = self.__dict_labels[parameter]
            except KeyError:
                label = "%s"
            labels.append(label)
        label = ", ".join(labels)
        label_size = len(label)
        if label_size < 2:
            label_font = 14
        elif label_size == 2:
            label_font = 12
        elif label_size > 2:
            label_font = 10

        # Checks that key is a valid measurement
        if key in self.measures is False:
            raise ValueError("%s is not in self.measures") % (key)
        else:
            pass

        # Tries to find sample name
        if "sample" in self.parameters:
            sample = self["sample"]
        else:
            pass

        # Looks for fig and ax
        try:
            fig = kwargs["fig"]
            ax = kwargs["ax"]
            kwargs.pop("fig")
            kwargs.pop("ax")
            return_fig = False
        except KeyError:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            return_fig = True

        zero_line = 0
        y_axis = None

        # Draws the curves
        if key == "Tp_Tm":
            p = [self[parameter] for parameter in parameters]
            ax.plot(self[x_axis], self["Tp"], label=r"T$^+$ "+label %
                    tuple(p), *args, **kwargs)
            ax.plot(self[x_axis], self["Tm"], label=r"T$^-$ "+label %
                    tuple(p), *args, **kwargs)
            ax.legend(fontsize=label_font)
        else:
            p = [self[parameter] for parameter in parameters]
            x_data = self[x_axis]
            if key in ["kxy", "kxy/T"]:
                y_data = 10*self[key]
            else:
                y_data = self[key]
            ax.plot(self[x_axis], self[key], label=label %
                    tuple(p), *args, **kwargs)
            if self[key].min()*self[key].max() < 0 and zero_line == 0:
                ax.plot(self[x_axis], 0*self[key], "--k", lw=2)
                zero_line += 1

            elif self[key].min()*self[key].max() > 1 and zero_line == 0:
                if self[key].max() < 0:
                    y_axis = "Negative"
                else:
                    y_axis = "Positive"

        # Makes it pretty
        ax.set_xlabel(self.__dict_axis[x_axis], fontsize=axis_fs)
        ax.set_ylabel(self.__dict_axis[key], fontsize=axis_fs)
        ax.tick_params(axis="both", which="both", direction="in",
                       top=True, right=True)

        if label_size != 0:
            ax.legend(fontsize=label_font)
        else:
            pass

        # If sample is the same for all measurements print it on the figure
        if len(fig.axes) == 1:
            fig.tight_layout(rect=[0.01, 0.01, 1, 0.95])
            plt.figtext(0.05, 0.005, sample, fontsize=axis_fs -
                        2, va="bottom", ha="left")
        else:
            pass

        # Set axis to start at 0
        if x_axis in ["T_av", "T0"]:
            ax.set_xlim(0)
        else:
            pass
        if y_axis is not None:
            if y_axis == "Positive":
                ax.set_ylim(0, ax.get_ylim()[1])
            else:
                ax.set_ylim(ax.get_ylim()[0], 0)
        else:
            pass
        if show is False:
            plt.close()
        elif show is True:
            plt.show()

        if return_fig is True:
            return fig, ax
        else:
            return

    def Plot_all(self, *args, **kwargs):
        """
        Plots all non trivial measures, all the same kwargs as Conductivity.Plot
        with the addition of filename to save the file.
        """

        remove = ["T_av", "T0", "Tp", "Tm", "T0_fit", "I_fit"]
        measures = [i for i in self.measures if i not in remove]
        figures = []

        try:
            save = kwargs["save"]
            kwargs.pop("save")
        except KeyError:
            save = False

        try:
            filename = kwargs["filename"]
            kwargs.pop("filename")
        except KeyError:
            if save is True:
                filename = self["filename"].replace(".dat", ".pdf")
                filename = filename.replace("data", "figures")
                directory = os.path.split(filename)[0]
                print(directory)

                # Creates a figure directory with the same structure as the data
                # directory if the answer is yes, otherwise saves the pdf file
                # with the data
                if os.path.isdir(directory) is False:
                    answer = input(
                        "Do you want to create directory (Y/n): %s" % directory)
                    if answer in ["Y", "y", "", "yes"]:
                        os.makedirs(directory)
                    else:
                        filename = self["filename"].replace(".dat", ".pdf")
                        print("Figures will be saved with data at %s" %
                              filename)
            else:
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

    def Plot_fancy(self, *args, **kwargs):
        """
        Just like Plot_all but with a fancy layout that is more suited to
        ipython notebooks
        """

        remove = ["T_av", "T0", "Tp", "Tm", "I_fit", "T0_fit"]

        measures = [i for i in self.measures if i not in remove]
        ref_meas = ["kxx", "kxx/T", "kxy", "kxy/kxx", "dTx",
                    "dTx/T", "dTy", "dTy/dTx", "Resistance", "Tp_Tm"]
        measures = [i for i in ref_meas if i in measures]

        fig, ax = self.__create_grid(measures)

        try:
            kwargs.pop("show")
        except KeyError:
            pass

        try:
            filename = kwargs["filename"]
            kwargs.pop("filename")
        except KeyError:
            filename = None

        for i in range(len(measures)):
            self.Plot(measures[i], *args, show=None,
                      fig=fig, ax=ax[i], **kwargs)

        if hasattr(self, "__sample") is True:
            plt.suptitle(self["sample"], fontsize=22)
        else:
            pass

        fig.tight_layout(rect=[0.01, 0.01, 1, 0.95])

        if filename is not None:
            filename = os.path.abspath(filename)
            pp = PdfPages(filename)
            pp.savefig(fig)
            pp.close()
        else:
            pass

        return fig, ax

    def Write_out(self, filename=None):
        """
        Writes the treated data to a file
        """
        if filename is None:
            if self["H"] == "0.0" or self["symmetrize"] is False:
                filename = self["filename"].replace(".dat", "-treated.dat")
            else:
                filename = self["filename"].replace(".dat", "-sym-treated.dat")
        else:
            filename = os.path.abspath(filename)

        parameters1 = ["sample", "date", "mount", "H"]
        parameters2 = ["w", "t", "L"]
        measures = ["T_av", "T0", "Tp", "Tm", "dTx", "kxx", "dTy", "kxy"]
        columns = ["T_av(K)", "T0(K)", "T+(K)", "T-(K)",
                   "dTx(K)", "kxx(W/Km)", "dTy(K)", "kxy(W/Km)"]
        if self["H"] == "0.0":
            measures = measures[0:6]
            columns = columns[0:6]
        else:
            pass

        columns = "\t".join(columns)

        comments1 = "\n".join(["%s\t=\t%s" % (i, self[i])
                               for i in parameters1])
        comments2 = "\n".join(["%s\t=\t%1.3e" % (i, self[i])
                               for i in parameters2])
        header = comments1+"\n"+comments2+"\n"+columns
        data = np.array([self[i] for i in measures]).T

        np.savetxt(filename, data, delimiter="\t",
                   header=header, fmt="%.6e")

    def Current(self, _min, _max, deg=5, T_max=100, N=100, *args, **kwargs):
        """
        Used to compute the optimal current function for the sample.

        Parameters:
        ------------------------------------------------------------------------
        _min, _max: int or float
                The min/max of dT/T in percentages

        deg:        int
                The degree of the polynomial fit

        T_max:      int or float
                T0 max for the plot

        N:          int
                Number of points in the plot
        """

        datafile = os.path.abspath("/".join(self["filename"].split("/")[0:-1]))
        rnge = "%1.0f%s_to_%1.0f%s.dat" % (_min, "%", _max, "%")
        name = "_".join([self["sample"].replace(" ", "_"), "dTovT", rnge])
        datafile = os.path.join(datafile, name)
        n = self["T_av"].shape[0]
        dT_T = np.linspace(_min/100, _max/100, n)
        alpha = self["w"]*self["t"]/self["L"]
        I = np.sqrt(self["kxx"]*alpha*self["T_av"]*dT_T/5000)
        coeff_I = np.polyfit(self["T0"], I, deg)
        poly_func = np.poly1d(coeff_I)
        T0 = np.linspace(0, T_max, N)
        I_fit = poly_func(T0)

        self["T0_fit"] = T0
        self["I_fit"] = I_fit*1000
        self["coeff_I"] = coeff_I
        self.measures += ["T0_fit", "I_fit"]

        # Looks for show as kwarg
        try:
            show = kwargs["show"]
            if type(show) is not bool:
                if show is not None:
                    raise TypeError("show must be of type bool or None")
                else:
                    kwargs.pop("show")
            else:
                kwargs.pop("show")
        except KeyError:
            show = True

        try:
            filename = kwargs["filename"]
            kwargs.pop("filename")
        except KeyError:
            filename = None

        label = r"$\Delta$ T / T from %1.2f%s to %1.2f%s" % (
            _min, "%", _max, "%")

        fig, ax = self.Plot("I_fit", x_axis="T0_fit", show=None, parameters=[])
        plt.figtext(1-0.005, 0.005, label, fontsize=14,
                    va="baseline", ha="right")

        if show is True:
            plt.show()
        elif show is False:
            plt.close()
        else:
            pass

        if filename is not None:
            filename = os.path.abspath(filename)
            pp = PdfPages(filename)
            pp.savefig(fig)
            pp.close()
        else:
            pass

        try:
            write = kwargs["write"]
            kwargs.pop("write")
        except KeyError:
            if os.path.isfile(datafile) is False:
                write = True
            else:
                answer = input(
                    "File %s already exists, overwrite? (y/N)" % datafile)
                if answer in ["Y", "y", "O", "o", "yes", "Yes", "oui", "Oui"]:
                    write = True
                    print("File overwritten")
                else:
                    write = False
                    print("File will not be saved")

        if write is True:
            degrees = np.array([i for i in range(coeff_I.shape[0])])
            data = np.array([degrees, coeff_I[::-1]]).T
            header = "Current function coefficients\norder\tcoeff"
            np.savetxt(datafile, data, delimiter="\t",
                       header=header, fmt=["%i", "%.18e"])
        else:
            pass

        return

    def __getitem__(self, key):
        if type(key) is str:
            return getattr(self, "__"+key)
        else:
            C = Conductivity()

            for i in self.raw_data:
                setattr(C, "__"+i, getattr(self, "__"+i)[key])

            for i in self.measures:
                if i != "Tp_Tm":
                    setattr(C, "__"+i, getattr(self, "__"+i)[key])
                else:
                    setattr(C, "__"+i, None)
            for i in self.parameters:
                setattr(C, "__"+i, getattr(self, "__"+i))

            misc = self.__dict__.keys()
            for k in misc:
                if hasattr(C, k) is True:
                    pass
                else:
                    setattr(C, k, getattr(self, k))

            C.measures = self.measures
            C.parameters = self.parameters
            return C

    def __setitem__(self, key, value):
        if type(key) is str:
            setattr(self, "__"+key, value)
        else:
            pass
        return

    def __delitem__(self, key):
        delattr(self, "__"+key)
        return

    def Get_known_measures(self):
        return list(self.__dict_measures.keys())

################################################################################
#                       ____   ____ ____  ___ ____ _____                       #
#                      / ___| / ___|  _ \|_ _|  _ \_   _|                      #
#                      \___ \| |   | |_) || || |_) || |                        #
#                       ___) | |___|  _ < | ||  __/ | |                        #
#                      |____/ \____|_| \_\___|_|    |_|                        #
################################################################################


if __name__ == "__main__":

    try:
        filename = sys.argv[1]
    except IndexError:
        raise Exception("First argument must be a filename")
    try:
        w = float(sys.argv[2])
        t = float(sys.argv[3])
        L = float(sys.argv[4])
        sign = int(sys.argv[5])
        sample = Conductivity(filename, w, t, L, sign)
    except IndexError:
        try:
            w = float(sys.argv[2])
            t = float(sys.argv[3])
            L = float(sys.argv[4])
            sample = Conductivity(filename, w, t, L, sign)
        except IndexError:
            sample = Conductivity(filename)
    sample.Plot_all(save=True)
    sample.Write_out()
