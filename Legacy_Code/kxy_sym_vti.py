import sys
import os
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from thermometry import seebeck_thermometry

################################################################################
#                        _                _           _                        #
#                       / \   _ __   __ _| |_   _ ___(_)___                    #
#                      / _ \ | '_ \ / _` | | | | / __| / __|                   #
#                     / ___ \| | | | (_| | | |_| \__ \ \__ \                   #
#                    /_/   \_\_| |_|\__,_|_|\__, |___/_|___/                   #
#                                           |___/                              #
################################################################################

# Uses the first command line argument to determine
# where is the file containing the data for the positive field
data_file_up = os.path.abspath(sys.argv[1])

# Extracting the H field value from the filename
H = (data_file_up.split("/"))[-1].split("-")[2]  # Gets the string
h = float(H.split("T")[0])  # Convert the string to a float

# Gets the sample name from the file's comments
with open(data_file_up) as f:
    Sample = f.readline().strip().split(" ")[-1]

# Formated sample name + H field value to write on graphs
H_Sample = "H=%s" % H+" "+Sample

# Looks for a second file if the H field is not 0
if h != 0:
    data_file_down = data_file_up.replace("TS-", "TS--")
else:
    data_file_down = data_file_up

# Reads both files to an array
data_up = np.genfromtxt(data_file_up, delimiter="\t").T
data_down = np.genfromtxt(data_file_down, delimiter="\t").T

# Symetrization of the data
data = np.zeros(data_up.shape)
sym = [0, 1, 2, 5, 6, 7, 8, 9, 10]  # Indexes of symetric data
anti = [3, 4]  # Indexes of antisymetric data
data[sym] = 0.5*(data_up[sym]+data_down[sym])
data[anti] = 0.5*(data_up[anti]-data_down[anti])
del data_up, data_down

# Checks if geometry is specified
try:
    w = float(sys.argv[2])
    t = float(sys.argv[3])
    L = float(sys.argv[4])
except IndexError:
    w = 1e-6
    t = 1e-6
    L = 1e-6
alpha = w*t/L
# String to be added in output file header
alpha_string = "w = %1.3e\nt = %1.3e\nL = %1.3e\n" % (w, t, L)

T0 = data[0]
dT_abs = data[0]*0
dTx = data[0]*0
dT_prec = dT_abs+1000
# Boucle pour converger vers T réel
while abs(dT_prec.sum()-dT_abs.sum()) > 1e-10:
    dT_prec = dT_abs
    S1 = seebeck_thermometry((dT_abs/2+T0))
    S2 = seebeck_thermometry(T0+dT_abs+dTx/2)
    dT_abs = (data[6]-data[5])/S1/1000
    dTx = (data[8]-data[7])/S2/1000
    T_moins = dT_abs+T0
    T_plus = T_moins+dTx
    T_avg = T_moins+dTx/2

# Calculates the dissipated power
Q = 5000*data[2]*data[2]

# Calculates kxx
kxx = Q/dTx/alpha

#dTy and kxy
if h != 0:
    S = seebeck_thermometry((T_avg))
    dTy = -(data[4]-data[3])/S/1000
    kxy = kxx*dTy/dTx*L/w

################################################################################
#                         ___        _               _                         #
#                        / _ \ _   _| |_ _ __  _   _| |_                       #
#                       | | | | | | | __| '_ \| | | | __|                      #
#                       | |_| | |_| | |_| |_) | |_| | |_                       #
#                        \___/ \__,_|\__| .__/ \__,_|\__|                      #
#                                       |_|                                    #
#                                                                              #
################################################################################

    # Output array
# This code formats the array to be saved and generates a header containing
# the geometry of the sample and column titles
if h != 0:
    data_out = np.array([T_avg, T0, T_plus, T_moins, dTx, kxx, dTy, kxy]).T
    columns = ["T_av(K)", "T0(K)", "T+(K)", "T-(K)", "dTx(K)",
               "kxx(W/Km)", "dTy(K)", "kxy(W/Km)"]
    header = "\t".join(columns)
    header = alpha_string+header
else:
    data_out = np.array([T_avg, T0, T_plus, T_moins, dTx, kxx]).T
    columns = ["T_av(K)", "T0(K)", "T+(K)", "T-(K)", "dTx(K)", "kxx(W/Km)"]
    header = "\t".join(columns)
    header = alpha_string+header

filename = data_file_up.split(".")
filename[-2] += "_treated"
filename = ".".join(filename)
np.savetxt(filename, data_out, delimiter='\t', header=header,fmt="%.6e")

################################################################################
#                                               _                              #
#                          __ _ _ __ __ _ _ __ | |__  ___                      #
#                         / _` | '__/ _` | '_ \| '_ \/ __|                     #
#                        | (_| | | | (_| | |_) | | | \__ \                     #
#                         \__, |_|  \__,_| .__/|_| |_|___/                     #
#                         |___/          |_|                                   #
#                                                                              #
################################################################################

# Pdf file
pdf_file = data_file_up.replace("/data/", "/figures/")
pdf_file = pdf_file.replace(".dat", ".pdf")
directory = "/".join(pdf_file.split("/")[0:-1])

# Creates a figure directory with the same structure as the data
# directory if the answer is yes, otherwise saves the pdf file with
# the data
if os.path.isdir(directory) is False:
    answer = input("Do you want to create directory (Y/n): %s" % directory)
    if answer in ["Y", "y", "", "yes"]:
        os.makedirs(directory)
    else:
        pdf_file = data_file_up.replace(".dat", ".pdf")

# Pdf object to store the figures
pp = PdfPages(pdf_file)

# rcParams to make sure graphs are beautiful
mp.rc("xtick", direction="in", top=True)
mp.rc("ytick", direction="in", right=True)
mp.rc("figure", figsize=(7, 4))

# If statement to compute only when usefull
if h != 0:
    # kxy
    fig, ax = plt.subplots()
    ax.plot(T_avg, kxy*1000, "--o", c="green")
    if kxy.min()*kxy.max() < 0:
        ax.axhline(0, ls="--", c="k", lw="0.75")
    ax.set_xlabel("T (K)", fontsize=16)
    ax.set_ylabel(r"$\kappa_{xy}$ (mW/Km)", fontsize=16)
    ax.set_xlim(0)
    plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
    pp.savefig()

    # dTy
    fig, ax = plt.subplots()
    ax.plot(T_avg, dTy*1000, "--o", c="green")
    if dTy.min()*dTy.max() < 0:
        ax.axhline(0, ls="--", c="k", lw="0.75")
    ax.set_xlabel("T (K)", fontsize=16)
    ax.set_ylabel(r"$\Delta T_y$ (mK)", fontsize=16)
    ax.set_xlim(0)
    plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
    pp.savefig()

# kxx
fig, ax = plt.subplots()
ax.plot(T_avg, kxx, "--o", c="blue")
ax.set_xlabel("T(K)", fontsize=16)
ax.set_ylabel(r"$\kappa_{xx}$ (W/Km)", fontsize=16)
ax.set_xlim(0)
plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
pp.savefig()

# dTx
fig, ax = plt.subplots()
ax.plot(T_avg, dTx, "--o", c="gray")
ax.set_xlabel("T (K)", fontsize=16)
ax.set_ylabel(r"$\Delta T_x$ (K)", fontsize=16)
ax.set_xlim(0)
plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
pp.savefig()

# If statement to compute only of usefull
if h != 0:
    # dTy/dTx
    fig, ax = plt.subplots()
    ax.plot(T_avg, 100*dTy/dTx, '--o', c="grey")
    if (dTy/dTx).min()*(dTy/dTx).max() < 0:
        ax.axhline(0, ls="--", c="k", lw="0.75")
    ax.set_xlabel("T (K)", fontsize=16)
    ax.set_ylabel(r"$\Delta T_y/\Delta T_{x}$ (%)", fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
    ax.set_xlim(0)
    pp.savefig()

    # kxy/kxx
    fig, ax = plt.subplots()
    ax.plot(T_avg, 100*kxy/kxx, '--o', c="grey")
    if (kxy/kxx).min()*(kxy/kxx).max() < 0:
        ax.axhline(0, ls="--", c="k", lw="0.75")
    ax.set_xlabel("T (K)", fontsize=16)
    ax.set_ylabel(r"$\kappa_{xy}/\kappa_{xx}$ (%)", fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
    ax.set_xlim(0)
    pp.savefig()

# T+ T-
fig, ax = plt.subplots()
ax.plot(T_avg, T_plus, "--o", c="red", label=r"$T_+$")
ax.plot(T_avg, T_moins, "--o", c="blue", label="$T_-$")
ax.legend(fontsize=16)
ax.set_xlabel("T (K)", fontsize=16)
ax.set_ylabel("T (K)", fontsize=16)
plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
ax.set_xlim(0)
pp.savefig()

# dTx/T_avg
fig, ax = plt.subplots()
ax.plot(T_avg, dTx/T_avg*100, '--o', c="grey")
ax.set_xlabel("T (K)", fontsize=16)
ax.set_ylabel(r"$\Delta T_x/T_{avg}$ (%)", fontsize=16)
ax.set_xlim(0)
plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
pp.savefig()

#T_av-T0 / dTx
fig, ax = plt.subplots()
ax.plot(T_avg, (T_avg-T0)/dTx, '--o', c="grey")
ax.axhline(0.5, ls="--", c="r", lw="0.75",
           label="Ideal thermal contact resistance 1/2")
ax.set_ylabel(r"($T_{avg}-T_0$)/$\Delta T_x$", fontsize=16)
ax.set_xlabel("T (K)", fontsize=16)
ax.legend(frameon=False)
ax.set_xlim(0)
plt.figtext(0.025, 0.025, H_Sample, fontsize=14)
pp.savefig()

pp.close()
