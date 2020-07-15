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
    __dict_parameters["mount"] = ["mount"]
    __dict_parameters["sample"] = ["Sample", "sample"]
    __dict_parameters["date"] = ["Date", "date"]

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

    # Same principle then before but for curve labels
    __dict_labels = dict()
    __dict_labels["H"] = r"H = %sT"
    __dict_labels["sample"] = r"Sample: %s"
    __dict_labels["date"] = r"%s"

    def __init__(self, filename=None, w=1e-6, t=1e-6, L=1e-6, sign=1, **kwargs):

        self.parameters = []
        dict_geo = {"w": w, "t": t, "L": L}
        for key, value in dict_geo.items():
            setattr(self, "__"+key, value)
            self.parameters.append(key)

        for key, value in kwargs.items():
            setattr(self, "__"+key, value)
            self.parameters.append(key)

        if sign in [1, -1]:
            self["sign"] = sign
        else:
            raise ValueError("Sign must be 1 or -1")

        self.__read_file(filename)

        if getattr(self, "__H") != "0.0":
            self.__symetrize()

        self.__analyze()
        self.__add_measure()

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
                self["dTy"] = self["sign"]*(self["dTy_Q"]-self["dTy_0"])/1000/S
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
                self["dTy"] = self["sign"]*(self["dTy_Q"]-self["dTy_0"])/S/1000
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
        self["filename"] = filename

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

        # If sample is the same for all measurements print it on the figure
        if len(fig.axes) == 1:
            fig.tight_layout()
            plt.figtext(0.05, 0.005, sample, fontsize=axis_fs -
                        2, va="baseline", ha="left")
        else:
            pass

        # Makes it pretty
        ax.set_xlabel(self.__dict_axis[x_axis], fontsize=axis_fs)
        ax.set_ylabel(self.__dict_axis[key], fontsize=axis_fs)
        ax.tick_params(axis="both", which="both", direction="in",
                       top=True, right=True)

        if label_size != 0:
            ax.legend(fontsize=label_font)
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

        return fig, ax

    def Plot_all(self, *args, **kwargs):
        """
        Plots all non trivial measures, all the same kwargs as Conductivity.Plot
        with the addition of filename to save the file.
        """

        remove = ["T_av", "T0", "Tp", "Tm"]
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
                directory = "/".join(filename.split("/")[0:-1])
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

        remove = ["T_av", "T0", "Tp", "Tm"]

        measures = [i for i in self.measures if i not in remove]
        ref_meas = ["kxx", "kxx/T", "kxy", "kxy/T", "dTx",
                    "dTx/T", "dTy/dTx", "Resistance","Tp_Tm"]
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

        if hasattr(self,"__sample") is True:
            plt.suptitle(self["sample"], fontsize=22)
        else:
            pass

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

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
            filename = self["filename"].replace(".dat", "_treated.dat")
        else:
            pass

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
