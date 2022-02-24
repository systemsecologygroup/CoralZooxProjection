#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 20:43:39 2022

@author: araharin
"""

# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load, std, log
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, floor, arange
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np


"File used to produce Gross growth vs Cost of symbiosis """
#### Plotting data collected from ode solver #####
font = {'family' : 'normal',
        'weight' : 'bold'}
#matplotlib.rc('font',**{'family':'Times', 'sans-serif':['Times']})
#matplotlib.rc('text', usetex = True)
#matplotlib.rc('font', **font)


plt.rcParams["font.family"] = "arial"   

font0 = FontProperties()
font = font0.copy()
weight = "bold"
font.set_weight(weight)   

"""Parameters and functions"""
#### Adaptive dynamics model for coral-algae symbiosis (with bleaching), Time is in month #####
G_C = 10/12  # this is G_max in the model
a = 1.0768/12 # Symbiont specific growth rate - linear growth rate
b = 0.0633 # Exponential Growth constant of the symbiont


alpha = 1e-3 # coral investment linear cost parameter 

M_C = (10/12)*1e-3 # coral mortality 

# Regional coral carrying capacity
K_C_GBR = 438718.902336e10 # in cm^2
K_C_SEA = 1191704.68901e10 # in cm^2
K_C_CAR = 190374.697987e10  # in cm^2

K_C_List = [K_C_GBR, K_C_SEA, K_C_CAR]


Ksmax = 3e6 # healthy measure of carying capacity of symbiont per host biomass Gamma

beta = (12)*1e2

gammaH = 1e6 # minimun symbiont density found on healthy coral colonie
fsize = 18

r = (12)*1e3
error = 1e-2 # prevent division by zeros these scales are alright because we are dealing with very high numbers
error2 = 1e-2 # prevent division by zeros
error3 = 1e-2 # prevent division by zeros

# Temperature Parameter resp. for the different region GBR, SEA , CAR
T0_list = array([26.78, 28.13, 27.10]) 
skew_list = array([0.0002, 3.82, 1.06])
rho_list = array([1.0, 0.81, 0.89])

def GrossG(coral, symbiont, trait, SST, T0, rho, skew, K_C_Reg):
    TempCORListCenter = (TempList - T0)/rho
    Tcenter = (SST-T0)/rho
    NormCor = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)
    Gx1Forcing = G_C*norm.pdf(Tcenter)*norm.cdf(Tcenter*skew)/max(NormCor) 
    E = 1 - exp(-beta*trait)
    symbiontH = gammaH*coral
    kappa = symbiont/(symbiontH + symbiont + error3)
    Benefit = Gx1Forcing*kappa*E*(1-coral/K_C_Reg)
    return Benefit

def SymbCost(coral, symbiont, trait):
    K_symb = Ksmax*coral 
    cost_gam = symbiont/(K_symb + error)
    Cost =alpha*exp(r*trait)*(cost_gam)
    return Cost



"""Plot results"""
##### Whole ranges for the parameter N (Used to get Figure 3 in manuscript) 

scale = 1e-11
rawNum_GBR = arange(0.01, 0.1, 0.00005)
rawNum_SEA = arange(0.01, 0.1, 0.00005)
rawNum_CAR = arange(0.01, 0.1, 0.00005) 

fsize = 14 #18 #22
fsize2 = 16 #18

# Open Files and plot

Locations = ["GBR", "SEA", "CAR"]
RCP = ["RCP26", "RCP45", "RCP85"] # for filename
RCP_title =  ["RCP 2.6", "RCP 4.5", "RCP 8.5"]
Color_list = [(0.306, 0.431+0.1, 0.545+0.3), (0.839+0.025, 0.60+0.05, 0.20), (0.839+0.16, 0.363, 0.35)]
Fig_lab = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
#Fig_lab = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]

# Colorbar setting
import matplotlib.colors as mcolors
import matplotlib as mpl
 
fig =plt.figure(figsize=(40, 40))  # (14, 20)
ax = plt.subplot(1, 1, 1)

    
AddTime = 2000*12  # spine up of 2000 years

ParamsPlus= ["Gmax-Plus", "a-Plus", "b-Plus", "KC-Plus", "Ksmax-Plus", "MC-Plus", "alpha-Plus", "r-Plus", "beta-Plus", "GammaH-Plus"]
ind = 1

"""
cmx = mpl.cm
# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmapC = mpl.cm.Vega20c_r#cm.copper_r
cNorm = mpl.colors.Normalize(vmin=N_List[3]-2.5*1e-10, vmax=N_List[12])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmapC)
"""

# Plot effect of Forcing, plot with cummulative number of bleaching

TempList = linspace(0, 75, 1500)

count = 0

# Regional coral carrying capacity
K_C_GBR = 438718.902336e10 # in cm^2
K_C_SEA = 1191704.68901e10 # in cm^2
K_C_CAR = 190374.697987e10  # in cm^2

Ksmax = 3e6 # healthy measure of carying capacity of symbiont per host biomass Gamma

file_list = ["Results_2/"]

#  Speed of adapations specific for each region # index in the whole N ranges
GBR_N_index = 908 # index in rawNum_GBR  0.0554
SEA_N_index = 330 # index in rawNum_SEA  0.0265
CAR_N_index = 275 # index in rawNum_CAR  0.02375
reg_N_index_true = [GBR_N_index, SEA_N_index, CAR_N_index]

#for the simulations with the right parameters (in folder Results/)
GBR_N_index = 1 #  0.0554
SEA_N_index = 1 #  0.0265
CAR_N_index = 1 #  0.02375
reg_N_index = [GBR_N_index, SEA_N_index, CAR_N_index]

for z in xrange(len(Locations)):
    sub1 = plt.subplot(4, 3, 1+count)
    plt.title(Locations_title[z], fontsize = fsize+2, fontproperties = font)
    part1 = sub1.twinx()
    for v in xrange(2,3):#len(RCP)):
        rcp = RCP[v]
        
        file0 = open("Monthly-SST-scenarios/Months-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        time0 = load(file0, allow_pickle = True)
        file0.close()
        time = concatenate((arange(min(time0)-(AddTime+12)/12, min(time0), 1/12), time0))   # time are already given in years with months in decimals    

        
        file1sst = open("Monthly-SST-scenarios/SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        SST = load(file1sst, allow_pickle = True)
        file1sst.close()
        
        file2 = open("Results/CORAL-"+rcp+"-"+Locations[z]+".dat", "r")
        HOSTSet1 = load(file2, allow_pickle = True)
        file2.close()
        
        file3 = open("Results/TRAIT-"+rcp+"-"+Locations[z]+".dat", "r")
        TRAITSet1 = load(file3, allow_pickle = True)
        file3.close()
        
        file4 =  open("Results/SYMB-"+rcp+"-"+Locations[z]+".dat", "r")
        SYMBSet1 = load(file4, allow_pickle = True)
        file4.close()
        
        HOST = HOSTSet1
        SYMB = SYMBSet1
        TRAIT = TRAITSet1
        
        if Locations[z] == "GBR":
            N_index = reg_N_index[z]
            N_index_true = reg_N_index_true[z]

        elif Locations[z] == "SEA":
            N_index = reg_N_index[z]                               
            N_index_true = reg_N_index_true[z]
    
        elif Locations[z] == "CAR":
            N_index = reg_N_index[z]
            N_index_true = reg_N_index_true[z]
             
        Host = HOST[N_index][list(time).index(min(time0)):]
        Trait = TRAIT[N_index][list(time).index(min(time0)):]
        Symb = SYMB[N_index][list(time).index(min(time0)):] 
        
        GG = GrossG(Host, Symb, Trait, SST, T0_list[z], rho_list[z], skew_list[z], K_C_List[z])
        SC = SymbCost(Host, Symb, Trait)
        Norm = max(max(GG), max(SC))
        sub1.plot(time0, GG/Norm, linewidth = 2, color = "blue", label = "Gross growth", alpha = 0.65)
        sub1.plot(time0, -SC/Norm, linewidth = 2, color = "orange", label = "Cost of symbiosis", alpha = 0.75)
        
        

      
    if z == 2:
        sub1.legend(fontsize = fsize2-2)
     
    
    tcks = np.arange(-0.5, 1.1+0.2, 0.4)

    if z == 0:
        sub1.set_yticks(tcks)
        sub1.set_yticklabels(["%.1f"%tc for tc in tcks], fontsize = fsize2)
        sub1.set_ylabel("Magnitude \n(normalized)", fontsize = fsize2)
    else:
        sub1.set_yticks(tcks)
        sub1.set_yticklabels(["" for tc in tcks], fontsize = fsize2)
    
    part1.set_yticks(tcks)
    part1.set_yticklabels(["" for tc in tcks], fontsize = fsize2)

    sub1.set_ylim((-0.5, 1.1))
    part1.set_ylim((-0.5, 1.1))
    sub1.set_xticks([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100])
    sub1.set_xticklabels([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100], rotation = 45, fontsize = fsize2)
    plt.xlim((2010, 2100))
    
    plt.xlabel("Years", fontsize = fsize2)
    count +=1
# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()


plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.125,
right=0.825,
hspace=0.2,
wspace=0.14)


plt.show()
