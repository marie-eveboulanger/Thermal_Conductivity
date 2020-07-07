import sys
import os
import re
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from thermometry import seebeck_thermometry

#Uses the first command line argument to determine 
#where is the file containing the data for the positive
#field and reading it to an array
data_file_up = sys.argv[1]
H = (data_file_up.split("/"))[-1].split("-")[2]#Gets the value of H
h = int(float(H.split("T")[0]))


#Looks for a second file if H is not 0
if h != 0:
    data_file_down = data_file_up.replace("TS-","TS--")
else:
    data_file_down = data_file_up
data_up = np.genfromtxt(data_file_up,delimiter="\t").T
data_down = np.genfromtxt(data_file_down,delimiter="\t").T

#Gets the sample name
with open(data_file_up) as f:
    Sample = f.readline().strip().split(" ")[-1]

H_Sample = "H=%s"%H+" "+Sample

#Checks if alpha is specified
try:
    alpha = sys.argv[2]
except IndexError:
    alpha = 1e-6

#Polynomial fit over Heat-off resistances
coeff_plus_up = poly.polyfit(np.log(data_up[3]),np.log(data_up[0]),8)
coeff_moins_up = poly.polyfit(np.log(data_up[5]),np.log(data_up[0]),8)

if h != 0:
    coeff_plus_down = poly.polyfit(np.log(data_down[3]),np.log(data_down[0]),8)
    coeff_moins_down = poly.polyfit(np.log(data_down[5]),np.log(data_down[0]),8)

#Determines which data points cannot be calibrated and removes them
index_plus = np.where(data_up[4]<data_up[3][-1])
data_up = np.delete(data_up,index_plus,axis=-1)

if h != 0:
    data_down = np.delete(data_down,index_plus,axis=-1)

#Calculates T+ T- 
T_plus_up = np.exp(poly.polyval(np.log(data_up[4]),coeff_plus_up))
T_moins_up = np.exp(poly.polyval(np.log(data_up[6]),coeff_moins_up))

if h != 0:
    T_plus_down = np.exp(poly.polyval(np.log(data_down[4]),coeff_plus_down))
    T_moins_down = np.exp(poly.polyval(np.log(data_down[6]),coeff_moins_down))

#Symmetrization 
if h !=0:
    T0 = (data_up[0]+data_down[0])/2
    T_plus = (T_plus_up+T_plus_down)/2
    T_moins = (T_moins_up+T_moins_down)/2
else:
    T0 = data_up[0]
    T_plus = T_plus_up
    T_moins = T_moins_up

dTx = T_plus-T_moins
T_avg = 0.5*(T_plus+T_moins)

#Calculates the dissipated power
Q = 5000*data_up[2]*data_up[2]

#Calculates kxx
kxx = Q/dTx/1e-6

#dTy and kxy
if h != 0:
    dTy_off = (data_up[7]-data_down[7])/2/seebeck_thermometry((T_avg+T0)/2)
    dTy = (data_up[8]-data_down[8])/2/seebeck_thermometry((T_avg+T0)/2)-dTy_off
    kxy = -kxx*dTy/dTx

#Output array
if h != 0:
    data_out = np.array([T_avg,T0,T_plus,T_moins,dTx,kxx,dTy,kxy]).T
    header = "\t".join(["T_av","T0","T+","T-","dTx","kxx","dTy","kxy"])
else:
    data_out = np.array([T_avg,T0,T_plus,T_moins,dTx,kxx]).T
    header = "\t".join(["T_av","T0","T+","T-","dTx","kxx"])

filename = data_file_up.split(".")
filename[-2] += "_treated"
filename = ".".join(filename)
np.savetxt(filename,data_out,delimiter='\t',header=header)

#Pdf file
pdf_file = data_file_up.replace("/data/","/figures/")
pdf_file = pdf_file.replace(".dat",".pdf")
directory = "/".join(pdf_file.split("/")[0:-1])
if os.path.isdir(directory) is False:
    os.makedirs(directory)
pp = PdfPages(pdf_file)

#Plots the data

#rcParams
mp.rc("xtick",direction="in",top=True)
mp.rc("ytick",direction="in",right=True)

if h != 0:
    #kxy
    fig,ax = plt.subplots()
    ax.plot(T0,kxy,"--o",c="green")
    ax.set_xlabel("T (K)",fontsize=16)
    ax.set_ylabel(r"$\kappa_{xy}$ (mW/Km)",fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
    pp.savefig()

    #dTy
    fig,ax = plt.subplots()
    ax.plot(T0,dTy,"--o",c="green")
    ax.set_xlabel("T (K)",fontsize=16)
    ax.set_ylabel(r"$\Delta T_y$ (mK)",fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
    pp.savefig()

#kxx
fig,ax = plt.subplots()
ax.plot(T0,kxx,"--o",c="blue")
ax.set_xlabel("T(K)",fontsize=16)
ax.set_ylabel(r"$\kappa_{xx}$ (W/Km)",fontsize=16)
plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
pp.savefig()

#dTx
fig,ax = plt.subplots()
ax.plot(T0,dTx,"--o",c="gray")
ax.set_xlabel("T (K)",fontsize=16)
ax.set_ylabel(r"$\Delta T_x$ (K)",fontsize=16)
plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
pp.savefig()

if h != 0:
    #dTy/dTx
    fig,ax = plt.subplots()
    ax.plot(T0,dTy/10/dTx,'--o',c="grey")
    ax.set_xlabel("T (K)",fontsize=16)
    ax.set_ylabel(r"$\Delta T_y/\Delta T_{x}$ (%)",fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
    pp.savefig()

    #kxy/kxx
    fig,ax = plt.subplots()
    ax.plot(T0,kxy/10/kxx,'--o',c="grey")
    ax.set_xlabel("T (K)",fontsize=16)
    ax.set_ylabel(r"$\kappa_{xy}/\kappa_{xx}$ (%)",fontsize=16)
    plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
    pp.savefig()

#T+ T-
fig,ax = plt.subplots()
ax.plot(T0,T_plus,"--o",c="red",label=r"$T_+$")
ax.plot(T0,T_moins,"--o",c="blue",label="$T_-$")
ax.legend(fontsize=16)
ax.set_xlabel("T (K)",fontsize=16)
ax.set_ylabel("T (K)",fontsize=16)
plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
pp.savefig()

#dTx/T_avg
fig,ax = plt.subplots()
ax.plot(T0,dTx/T_avg*100,'--o',c="grey")
ax.set_xlabel("T (K)",fontsize=16)
ax.set_ylabel(r"$\Delta T_x/T_{avg}$ (%)",fontsize=16)
plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
pp.savefig()

#T_av-T0 / dTx
fig,ax = plt.subplots()
ax.plot(T0,(T_avg-T0)/dTx,'--o',c="grey")
ax.set_ylabel(r"($T_{avg}-T_0$)/$\Delta T_x$",fontsize=16)
ax.set_xlabel("T (K)",fontsize=16)
plt.figtext(0.025, 0.025, H_Sample,fontsize=14)
pp.savefig()

if h != 0:
    #Cernox
    fig,ax = plt.subplots(2,1)
    ax[0].plot(np.log(data_up[3]),poly.polyval(np.log(data_up[3]),coeff_plus_up),'-r',label=r"Fit $T_+$ $\dot{Q}$=0")
    ax[0].plot(np.log(data_up[4]),np.log(T_plus_up),'or',label=r"$T_+$ $\dot{Q}\neq0$")
    ax[0].plot(np.log(data_up[5]),poly.polyval(np.log(data_up[5]),coeff_moins_up),'-b',label=r"Fit $T_+$ $\dot{Q}$=0")
    ax[0].plot(np.log(data_up[6]),np.log(T_moins_up),'ob',label=r"$T_-$ $\dot{Q}\neq0$")
    ax[0].set_xlabel(r"$\log(V)$")
    ax[0].set_ylabel(r"$\log(T_0)$")
    ax[0].legend()
    ax[0].set_title("H=%s"%H)


    ax[1].plot(np.log(data_down[3]),poly.polyval(np.log(data_down[3]),coeff_plus_down),'-r',label=r"Fit $T_+$ $\dot{Q}$=0")
    ax[1].plot(np.log(data_down[4]),np.log(T_plus_down),'or',label=r"$T_+$ $\dot{Q}\neq0$")
    ax[1].plot(np.log(data_down[5]),poly.polyval(np.log(data_down[5]),coeff_moins_down),'-b',label=r"Fit $T_+$ $\dot{Q}$=0")
    ax[1].plot(np.log(data_down[6]),np.log(T_moins_down),'ob',label=r"$T_-$ $\dot{Q}\neq0$")
    ax[1].set_xlabel(r"$\log(V)$")
    ax[1].set_ylabel(r"$\log(T_0)$")
    ax[1].legend()
    ax[1].set_title("H=-%s"%H)
    plt.figtext(0.05, 0.5, H_Sample,fontsize=14)

    fig.tight_layout()

    pp.savefig()

pp.close()
