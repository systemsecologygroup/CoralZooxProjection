# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load, std, log
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, floor, arange
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np


"""File used to produce Figure 3 of manuscript"""

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

# Parameter resp. for the different region GBR, SEA , CAR
T0_list = array([26.78, 28.13, 27.10]) 
skew_list = array([0.0002, 3.82, 1.06])
rho_list = array([1.0, 0.81, 0.89])

# For model with bleaching
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

startTime = 2010 #1970

# Regional coral carrying capacity
K_C_GBR = 438718.902336e10 # in cm^2
K_C_SEA = 1191704.68901e10 # in cm^2
K_C_CAR = 190374.697987e10  # in cm^2

Ksmax = 3e6 # healthy measure of carying capacity of symbiont per host biomass Gamma

#  Speed of adapations specific for each region # index in the whole N ranges
GBR_N_index = 908 # index in rawNum_GBR  0.0554
SEA_N_index = 330 # index in rawNum_SEA  0.0265
CAR_N_index = 275 # index in rawNum_CAR  0.02375
reg_N_index_true = [GBR_N_index, SEA_N_index, CAR_N_index]

#for the simulations with the right parameters (in folder Results/)
GBR_N_index = 0 #  0.0554
SEA_N_index = 0 #  0.0265
CAR_N_index = 0 #  0.02375
reg_N_index = [GBR_N_index, SEA_N_index, CAR_N_index]

import pdb


for z in xrange(1):
    for v in xrange(1):
        rcp = RCP[v]
    
        file0 = open("Monthly-SST-scenarios/Months-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        time0 = load(file0, allow_pickle = True)
        file0.close()
        
        time = concatenate((arange(min(time0)-(AddTime+12)/12, min(time0), 1/12), time0))   # time are already given in years with months in decimals

        file1sst = open("Monthly-SST-scenarios/SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        SST = load(file1sst, allow_pickle = True)
        file1sst.close()
        
        file1 = open("Monthly-SST-scenarios/SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        TempNS = load(file1, allow_pickle = True)
        file1.close()
        
        sub1 = plt.subplot(3, 1, 1+count)
        #TraitTicks = array([-0.5, 0, 101-100, 102-100, 102.5-100])
        TraitTicks = array([-0.5, 0, 0.5, 1, 1.5, 2, 2.5])
        plt.xticks([time0[0], 2100], [time0[0], 2100], fontsize = 10)
        sub1.set_xlim((startTime, max(time0))) 
        part1 = sub1.twinx()     
        sub4 = plt.subplot(3, 1, 2+count)
          
        plt.xticks([time0[0], 2100], [time0[0], 2100], fontsize = 10)
        part4 = sub4.twinx()
        #hostTicks = array([0, 20, 40, 60, 80, 100, 120]) - 100    
        hostTicks = array([0-100, 20-100, 40-100, 60-100, 80-100, 0, 120-100])  
        
        sub7 = plt.subplot(3, 1, 3+count)
        plt.xticks([time0[0], 2100], [time0[0], 2100], fontsize = 10)

        Symbticks = array([0, 30, 60, 90, 120, 150, 180, 200]) - 100
        part7 = sub7.twinx()
        if count==3:
            sub7.set_ylabel("Symbiont abundance \n (% change)\n", fontsize = fsize + 2)
            sub7.text(startTime - 16, 180-100, Fig_lab[6+count], fontproperties=font, fontsize = fsize) # forcing
        else:
            sub7.text(startTime - 5, 180-100, Fig_lab[6+count], fontproperties=font, fontsize = fsize) # forcing
                    
        
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
        
            
        MeanlossCoralBiomass = []
        
        # in case one wants something specific from the runs in each regions, uncomment and modify
            
        if Locations[z] == "GBR":
            N_index = reg_N_index[z]
            N_index_true = reg_N_index_true[z]

        elif Locations[z] == "SEA":
            N_index = reg_N_index[z]                               
            N_index_true = reg_N_index_true[z]
    
        elif Locations[z] == "CAR":
            N_index = reg_N_index[z]
            N_index_true = reg_N_index_true[z]
             
        Host = HOST[N_index]
        Trait = TRAIT[N_index]
        Symb = SYMB[N_index]
        
        sub1.plot(time, Trait, linewidth = 2, color = "black", label = "eps = 0 (used)") 
        sub4.plot(time, Host, linewidth = 2, color = "black")
        sub7.plot(time, Symb, linewidth = 2, color = "black")
    
        for k in ("-0.2", "-0.1","+0.1", "+0.2"):
            file2 = open("Results_2/eps%s/"%k+"CORAL-"+rcp+"-"+Locations[z]+".dat", "r")
            HOSTSet1 = load(file2, allow_pickle = True)
            file2.close()
            
            file3 = open("Results_2/eps%s/"%k+"TRAIT-"+rcp+"-"+Locations[z]+".dat", "r")
            TRAITSet1 = load(file3, allow_pickle = True)
            file3.close()
            
            file4 =  open("Results_2/eps%s/"%k+"SYMB-"+rcp+"-"+Locations[z]+".dat", "r")
            SYMBSet1 = load(file4, allow_pickle = True)
            file4.close()
            
            HOST = HOSTSet1
            SYMB = SYMBSet1
            TRAIT = TRAITSet1
             
            Host = HOST[0]
            Trait = TRAIT[0]
            Symb = SYMB[0]
            
            sub1.plot(time, Trait, linewidth = 1, label = "eps = %s"%k) 
            sub4.plot(time, Host, linewidth = 1)
            sub7.plot(time, Symb, linewidth = 1)
    
        
        sub1.set_ylim((0, 5*np.std(Trait)))
        sub1.set_xlim((time[0], max(time0))) 
        part1.set_ylim((0, 5*np.std(Trait)))
        part1.set_xlim((time[0], max(time0))) 
        
        sub4.set_xlim((time[0], max(time0))) 
        sub4.set_ylim((0, 3.5*np.std(Host)))
        part4.set_ylim((0, 3.5*np.std(Host)))
        part4.set_xlim((time[0], max(time0))) 
        
        sub7.set_ylim((0, 3.5*np.std(Symb)))
        part7.set_ylim((0, 3.5*np.std(Symb))) 
        sub7.set_xlim((time[0], max(time0)))
        part7.set_xlim((time[0], max(time0))) 

        sub1.set_ylabel("Trait", fontsize = 18)
        sub4.set_ylabel("Coral abundance", fontsize = 18)
        sub7.set_ylabel("Symbiont abundance", fontsize = 18)
        
        #sub7.set_xlabel("Years", fontsize=fsize + 2)  
        sub7.text(time[100], -2e21, "The model dynamics before year %d is a spine-up phase \n which does not have any biological meaning \n but is an artefact an dynamical systems regardless of initial conditions"%time0[0], fontsize = 18)
        
        sub1.legend()
    count +=1 


# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

"""
plt.subplots_adjust(top=0.965,
bottom=0.095,
left=0.080,
#right=0.81,
right = 0.75,
hspace=0.175,
wspace=0.175)
"""
#plt.savefig("Fig3.pdf", bbox_inches = 'tight')    


plt.show()
