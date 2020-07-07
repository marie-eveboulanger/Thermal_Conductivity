import sys
import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from thermometry import seebeck_thermometry

#Uses the first command line argument to determine 
#where is the file containing the data and reading it to an array
data_file = sys.argv[1]
data = np.genfromtxt(data_file,delimiter="\t",skip_header=0).T

try:
    alpha = sys.argv[2]
except IndexError:
    alpha = 1e-6

#Polynomial fit over Heat-off resistances
coeff_plus = np.polyfit(data[3],data[0],15)
coeff_moins = np.polyfit(data[5],data[0],15)

#Determines which data points cannot be calibrated and removes them
index_plus = np.where(data[4]<data[3][-1])
data = np.delete(data,index_plus,axis=-1)

#Calculates various temperatures
T_plus = np.poly1d(coeff_plus)(data[4])
T_moins = np.poly1d(coeff_moins)(data[6])
T0 = data[0]
Delta_T = T_plus-T_moins
T_avg = 0.5*(T_plus+T_moins)

#Calculates the dissipated power
Q = 5000*data[2]*data[2]

#Calculates kxx
kxx = Q/Delta_T/1e-6

#Output array
data_out = np.array([T_avg,T0,T_plus,T_moins,Delta_T,kxx]).T
filename = data_file.split(".")
filename[-2] += "_treated"
filename = ".".join(filename)
header = "\t".join(["T_av","T0","T+","T-","dTx","kxx"])
np.savetxt(filename,data_out,delimiter='\t',header=header)

#Pdf file
pdf_file = data_file.replace("/data/","/figures/")
pdf_file = pdf_file.replace(".dat",".pdf")
directory = "/".join(pdf_file.split("/")[0:-1])
if os.path.isdir(directory) is False:
    os.makedirs(directory)
pp = PdfPages(pdf_file)

#Plots the data

#rcParams
mp.rc("xtick",direction="in",top=True)
mp.rc("ytick",direction="in",right=True)

#kxx
fig,ax = plt.subplots()
ax.plot(T0,kxx,"--o",c="blue")
ax.set_xlabel("T(K)",fontsize=16)
ax.set_ylabel(r"$\kappa_{xx}$ (W/Km)",fontsize=16)
pp.savefig()

#dTx
fig,ax = plt.subplots()
ax.plot(T0,Delta_T,"--o",c="gray")
ax.set_xlabel("T(K)",fontsize=16)
ax.set_ylabel(r"$\Delta T_x$ (K)",fontsize=16)
pp.savefig()

#T+ T-
fig,ax = plt.subplots()
ax.plot(T0,T_plus,"--o",c="red",label=r"$T_+$")
ax.plot(T0,T_moins,"--o",c="blue",label="$T_-$")
ax.legend(fontsize=16)
ax.set_xlabel("T (K)",fontsize=16)
ax.set_ylabel("T (K)",fontsize=16)
pp.savefig()

#dTx/T0
fig,ax = plt.subplots()
ax.plot(T0,Delta_T/T0*100,'--o',c="grey")
ax.set_xlabel("T (K)",fontsize=16)
ax.set_ylabel(r"$\Delta T_x/T_{0}$ (%)",fontsize=16)
pp.savefig()

#T_av-T0 / dTx
fig,ax = plt.subplots()
ax.plot(T_avg,(T_avg-T0)/Delta_T,'--o',c="grey")
ax.set_ylabel(r"($T_{avg}-T_0$)/$\Delta T_x$",fontsize=16)
ax.set_xlabel("T(K)",fontsize=16)
pp.savefig()

pp.close()
