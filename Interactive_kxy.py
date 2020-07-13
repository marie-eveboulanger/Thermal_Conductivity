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
    __dict_raw["T0"] = ["T0(K)"]
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
    __dict_parameters["Loc"] = ["BOT", "TOP", "Bot", "Top", "bot", "top"]
    __dict_parameters["sample"] = ["Sample", "sample"]
    __dict_parameters["date"] = ["Date", "date"]

    def __init__(self, data_file, w=1e-6, t=1e-6, L=1e-6, sign=1, raw_data=True):
        """
        Used to initialize the object

        Parameters:
        -----------------------------
        data_file:  string
        The location of the file
        Ex: /data/data_file.dat
        Only one file has to be specified even for measurements
        that include a magnetic field the code will find
        the matching file at the inverse field on its own.

        w:          float
        t:          float
        L:          float

        w,t,L are the dimensions of
        the sample

        raw_data:   boolean
        If raw_data is False the data_file should be
        a file containing the analyzed data that is outputed
        by the Write_data function or has the same structure.
        """
        data_file = os.path.abspath(data_file)
        self.sign = sign
        if raw_data == False:
            self.w = 1e-6
            self.t = 1e-6
            self.L = 1e-6
            self.data_out = np.genfromtxt(data_file, delimiter="\t").T
            with open(data_file) as df:
                for i, line in enumerate(df):
                    if i > 5:
                        break
                    else:
                        j = line.strip().split(" ")[0].split("\t")
                        if j[1] == "w":
                            self.w = float(j[-1])
                        elif j[1] == "t":
                            self.t = float(j[-1])
                        elif j[1] == "L":
                            self.L = float(j[-1])
                        elif j[1] == "H":
                            self.H = float(j[-1])
                            self.H_str = "%.1fT" % self.H
                        elif j[1] == "Sample":
                            self.sample = j[-1]
            self.data = None
            self.data_file = None
            return

        self.data = None
        self.data_out = None
        f = data_file.split("/")[-1]

        # If -H file is given
        if len(f.split("--")) == 2:
            # Finding the matching file at +H
            self.data_file_down = data_file
            self.data_file_up = data_file.replace("--", "-")
            self.data_file = None

            # Extracting the H field value from the filename
            self.H_str = (self.data_file_up.split("/"))[-1].split("-")[2]
            self.H = float(self.H_str.split("T")[0])
        # If +H file is given
        else:
            # Extracting the H field value from the filename
            self.H_str = (data_file.split("/"))[-1].split("-")[2]
            self.H = float(self.H_str.split("T")[0])
            if self.H == 0:
                self.data_file_up = None
                self.data_file_down = None
                self.data_file = data_file
            else:
                self.data_file_up = data_file
                self.data_file_down = data_file.replace("TS-", "TS--")
                self.data_file = None

        # Geometry
        self.w = w
        self.t = t
        self.L = L

        # Sample name and probe detection
        lines = []
        with open(data_file) as df:
            for i, line in enumerate(df):
                if i > 2:
                    break
                else:
                    lines.append(line)
        self.sample = lines[0].strip().split(" ")[-1]
        self.comments = lines[1]
        if lines[2].split("\t")[3] == "R+_0(V)":
            self.probe = "Tallahassee"
        else:
            self.probe = "VTI"
        return

    def __repr__(self):
        string = "Data for %s at H=%s" % (self.sample, self.H_str)
        return string

    def __call__(self):
        return

    def __read_file(self, filename):
        """
        Used to read the file header and the data. Also detects the probe that
        has been used for the measurement. Also detects if the data is raw or
        already analyzed.
        """

        # Converting filename to an absolute path if it is relative
        filename = os.path.abspath(filename)

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
        for line in self.lines:
            l = line.split("\t")
            # Checks for comments shorter than the raw data header
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
                for key, values in self.__dict_raw.items():
                    for i in range(len(l)):
                        if l[i].strip() in values:
                            setattr(self, "__"+key, raw_data[i])
                            self.raw_data.append(key)
                        else:
                            pass
                        if len(self.measures) == 0:
                            check_treated = True
                        else:
                            check_treated = False

                if check_treated is True:
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


    def Symmetrize(self):
        """
        This function symmetrizes the data using data at
        both H and -H
        """
        if self.data_file is not None:
            raise Exception("H=0 no need to symmetrize the data")
        else:
            data_up = np.genfromtxt(self.data_file_up, delimiter="\t").T
            data_down = np.genfromtxt(self.data_file_down, delimiter="\t").T
            data = np.zeros(data_up.shape)

            if self.probe == "VTI":
                sym = [0, 1, 2, 5, 6, 7, 8, 9, 10]
                anti_sym = [3, 4]
            elif self.probe == "Tallahassee":
                sym = [0, 1, 2, 3, 4, 5, 6, 9, 10]
                anti_sym = [7, 8]

            data[sym] = 0.5*(data_up[sym]+data_down[sym])
            data[anti_sym] = 0.5*(data_up[anti_sym]-data_down[anti_sym])
            self.data = data
            return

    def Analyze(self):
        """
        This function annalyzes the data to return usefull
        info like T_av and kxx etc
        """
        if self.data is None:
            if self.H == 0:
                self.data = np.genfromtxt(self.data_file, delimiter="\t").T
            else:
                print("Symetrizing the data")
                self.Symmetrize()

        alpha = self.w*self.t/self.L

        if self.probe == "VTI":
            T0 = self.data[0]
            Q = 5000*self.data[2]**2
            dT_abs = self.data[0]*0
            dTx = self.data[0]*0
            T_av = dT_abs*0
            prev = dT_abs+1000
            # Loop to converge towards actual temperatures
            while abs(prev.sum()-T_av.sum()) > 1e-10:
                prev = T_av
                S1 = seebeck_thermometry((dT_abs/2+T0))
                S2 = seebeck_thermometry(T0+dT_abs+dTx/2)
                dT_abs = abs(self.data[6]-self.data[5])/S1/1000
                dTx = (self.data[8]-self.data[7])/S2/1000
                T_minus = dT_abs+T0
                T_plus = T_minus+dTx
                T_av = T_minus+dTx/2
            S = seebeck_thermometry(T_av)
            dTy = self.sign*(self.data[4]-self.data[3])/S/1000

        elif self.probe == "Tallahassee":
            # Polynomial fit over Heat-off resistances
            coeff_plus = npp.polyfit(
                np.log(self.data[3]), np.log(self.data[0]), 8)
            coeff_minus = npp.polyfit(
                np.log(self.data[5]), np.log(self.data[0]), 8)

            # Determines which data points cannot be calibrated
            index = np.where(self.data[4] < self.data[3][-1])
            data = np.delete(self.data, index, axis=-1)

            # Calculates Everything
            T0 = data[0]
            Q = 5000*data[2]**2
            T_plus = np.exp(npp.polyval(np.log(data[4]), coeff_plus))
            T_minus = np.exp(npp.polyval(np.log(data[6]), coeff_minus))
            dTx = T_plus-T_minus
            T_av = 0.5*(T_plus+T_minus)
            S = seebeck_thermometry((T_av+T0)/2)
            dTy = self.sign*(data[8]-data[7])/1000/S

        kxx = Q/dTx/alpha
        kxy = kxx*dTy/dTx*self.L/self.w

        data_out = np.array([T_av, T0, T_plus, T_minus, dTx, kxx, dTy, kxy])
        self.data_out = data_out
        return

    def Write_data(self, filename=None):
        """
        Writes the analyzed data to a file
        """
        if filename is None:
            if self.data_file is None:
                filename = self.data_file_up.split(".")
                filename[-2] += "_treated"
                filename = ".".join(filename)
            else:
                filename = self.data_file.split(".")
                filename[-2] += "_treated"
                filename = ".".join(filename)

        w, t, L = self.w, self.t, self.L
        alpha_str = "w\t=\t%1.3e\nt\t=\t%1.3e\nL\t=\t%1.3e\n" % (w, t, L)
        sample_str = "Sample\t:\t%s\n" % self.sample
        H_str = "H\t=\t%.2f\n" % self.H
        columns = ["T_av(K)", "T0(K)", "T+(K)", "T-(K)",
                   "dTx(K)", "kxx(W/Km)", "dTy(K)", "kxy(W/Km)"]

        if self.H == 0:
            data_out = self.data_out[0:6].T
            header = sample_str+H_str+alpha_str+"\t".join(columns[0:6])
        else:
            data_out = self.data_out.T
            header = sample_str+H_str+alpha_str+"\t".join(columns)

        np.savetxt(filename, data_out, delimiter="\t",
                   header=header, fmt="%.6e")
        return

    def Kxy(self, show=True):
        """
        Plots kxy versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, kxy = self.data_out[0], self.data_out[7]

        fig, ax = plt.subplots()
        ax.plot(T_av, kxy*1000, "--o", c="green")

        if kxy.min()*kxy.max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\kappa_{xy}$ (mW/Km)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        else:
            plt.close()
        return fig, ax

    def DTy(self, show=True):
        """
        Plots dTy versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, dTy = self.data_out[0], self.data_out[6]

        fig, ax = plt.subplots()
        ax.plot(T_av, dTy*1000, "--o", c="green")

        if dTy.min()*dTy.max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\Delta T_y$ (mK)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        else:
            plt.close()
        return fig, ax

    def Kxx(self, show=None):
        """
        Plots kxx versus T_avg
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, kxx = self.data_out[0], self.data_out[5]

        fig, ax = plt.subplots()
        ax.plot(T_av, kxx, "--o", c="blue")

        if kxx.min()*kxx.max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T(K)", fontsize=16)
        ax.set_ylabel(r"$\kappa_{xx}$ (W/Km)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def Kxx_over_T(self, show=None):
        """
        Plots kxx versus T_avg
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, kxx = self.data_out[0], self.data_out[5]

        fig, ax = plt.subplots()
        ax.plot(T_av, 10*kxx/T_av, "--o", c="blue")

        if kxx.min()*kxx.max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T(K)", fontsize=16)
        ax.set_ylabel(r"$\kappa_{xx}/T$ (mW/K$^2$cm)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def DTx(self, show=None):
        """
        Plots dTx versus t_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, dTx = self.data_out[0], self.data_out[4]

        fig, ax = plt.subplots()
        ax.plot(T_av, dTx, "--o", c="gray")

        if dTx.min()*dTx.max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\Delta T_x$ (K)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def DTy_over_DTx(self, show=None):
        """
        Plots Dty/Dtx*100 versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, dTx, dTy = self.data_out[0], self.data_out[4], self.data_out[6]

        fig, ax = plt.subplots()
        ax.plot(T_av, 100*dTy/dTx, '--o', c="grey")

        if (dTy/dTx).min()*(dTy/dTx).max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\Delta T_y/\Delta T_{x}$ (%)", fontsize=16)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)
        ax.set_xlim(0)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def Kxy_over_Kxx(self, show=None):
        """
        Plots Dty/Dtx*100 versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, kxx, kxy = self.data_out[0], self.data_out[5], self.data_out[7]

        fig, ax = plt.subplots()
        ax.plot(T_av, 100*kxy/kxx, '--o', c="grey")

        if (kxy/kxx).min()*(kxy/kxx).max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\kappa_{xy}/\kappa_{xx}$ (%)", fontsize=16)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)
        ax.set_xlim(0)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def Temperatures(self, show=None):
        """
        Plots T+ and T- versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av = self.data_out[0]
        T_plus, T_minus = self.data_out[2], self.data_out[3]

        fig, ax = plt.subplots()
        ax.plot(T_av, T_plus, "--o", c="red", label=r"$T_+$")
        ax.plot(T_av, T_minus, "--o", c="blue", label="$T_-$")
        ax.legend(fontsize=16)
        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel("T (K)", fontsize=16)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)
        ax.set_xlim(0)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def DTx_over_T_av(self, show=None):
        """
        Plots dTx/T_av*100 versus T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, dTx = self.data_out[0], self.data_out[4]

        fig, ax = plt.subplots()
        ax.plot(T_av, dTx/T_av*100, '--o', c="grey")

        if (dTx/T_av).min()*(dTx/T_av).max() < 0:
            ax.axhline(0, ls="--", c="k", lw="0.75")

        ax.set_xlabel("T (K)", fontsize=16)
        ax.set_ylabel(r"$\Delta T_x/T_{avg}$ (%)", fontsize=16)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def Thermal_contact_resistance(self, show=None):
        """
        Plots the thermal contact resistance vs T_av
        """
        if self.data_out is None:
            self.Analyze()

        # rcParams to make sure graphs are beautiful
        mp.rc("xtick", direction="in", top=True)
        mp.rc("ytick", direction="in", right=True)
        mp.rc("figure", figsize=(7, 4))

        T_av, T0, dTx = self.data_out[0], self.data_out[1], self.data_out[4]

        fig, ax = plt.subplots()
        ax.plot(T_av, abs(T_av-T0)/dTx, '--o', c="grey")
        ax.axhline(0.5, ls="--", c="r", lw="0.75",
                   label="Ideal thermal contact resistance 1/2")
        ax.set_ylabel(r"($T_{avg}-T_0$)/$\Delta T_x$", fontsize=16)
        ax.set_xlabel("T (K)", fontsize=16)
        ax.legend(frameon=False)
        ax.set_xlim(0)
        plt.figtext(0.1, 0.025, self.sample, fontsize=14)
        plt.figtext(0.8, 0.025, "H="+self.H_str, fontsize=14)

        if show == True:
            plt.show()
        elif show == False:
            plt.close()
        return fig, ax

    def Plot_all(self, filename=None, save=True, show=False):
        """
        Plots all availlable figures and saves them in
        a single pdf using matplotlib's PDF_Pages
        Parameters:
        ----------------------------------------------
        filename:   string
            Output file name and path default is the
            same as data_file or data_file_up except
            the directory /data/ will be swapped for
            /figures/

        show: Boolean
            Determines if the figures are shown on screen
        """
        if save == True:
            if filename is None:
                # Pdf file
                if self.data_file is None:
                    filename = self.data_file_up.replace("/data/", "/figures/")
                else:
                    filename = self.data_file.replace("/data/", "/figures/")
                filename = filename.replace(".dat", ".pdf")
                directory = "/".join(filename.split("/")[0:-1])

                # Creates a figure directory with the same structure as the data
                # directory if the answer is yes, otherwise saves the pdf file
                # with the data
                if os.path.isdir(directory) is False:
                    answer = input(
                        "Do you want to create directory (Y/n): %s" % directory)
                    if answer in ["Y", "y", "", "yes"]:
                        os.makedirs(directory)
                    else:
                        if self.data_file is None:
                            filename = self.data_file_up.replace(
                                ".dat", ".pdf")
                        else:
                            filename = self.data_file.replace(".dat", ".pdf")
                            print("Figures will be saved with data at %s" %
                                  filename)

            pp = PdfPages(filename)

        figures = []
        if self.H != 0:
            figures.append(self.Kxy(show)[0])
            figures.append(self.DTy(show)[0])
        figures.append(self.Kxx(show)[0])
        figures.append(self.Kxx_over_T(show)[0])
        figures.append(self.DTx(show)[0])
        if self.H != 0:
            figures.append(self.DTy_over_DTx(show)[0])
            figures.append(self.Kxy_over_Kxx(show)[0])
        figures.append(self.Temperatures(show)[0])
        figures.append(self.DTx_over_T_av(show)[0])
        figures.append(self.Thermal_contact_resistance(show)[0])

        if save is True:
            for i in figures:
                pp.savefig(i)
            pp.close()

        return figures

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
        sample = Conductivity(filename, w, t, L)
    except IndexError:
        sample = Conductivity(filename)
    sample.Plot_all()
    sample.Write_data()
