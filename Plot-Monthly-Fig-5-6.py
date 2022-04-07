# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load, std, log
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, floor, arange
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import pdb
plt.switch_backend('agg')

#### Plotting data collected from ode solver #####
font = {'family' : 'normal',
        'weight' : 'bold'}
#matplotlib.rc('font',**{'family':'Times', 'sans-serif':['Times']})
#matplotlib.rc('text', usetex = True)
#matplotlib.rc('font', **font)


"""
Sensitivity with respect to speed of acclimation Fig 5 and 6: Requires file generated Sensitivity-RCP-monthly-2.py"
"""

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
scale = 1e-11
rawNum_GBR = arange(0.01, 0.1, 0.00005)
rawNum_SEA = arange(0.01, 0.1, 0.00005)
rawNum_CAR = arange(0.01, 0.1, 0.00005) 

fsize = 12*4 #18 #22
fsize2 = 11*4 + 5 #18

# Open Files and plot

Locations = ["GBR", "SEA", "CAR"]
RCP = ["RCP26", "RCP45", "RCP85"] # for filename
RCP_title =  ["RCP 2.6", "RCP 4.5", "RCP 8.5"]
Color_list = [(0.651, 0.325, 0.529), (0.839+0.025, 0.60+0.05, 0.20), (0.839+0.16, 0.363, 0.35)]
Fig_lab = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
#Fig_lab = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]

# Colorbar setting
import matplotlib.colors as mcolors
import matplotlib as mpl
 
fig =plt.figure(figsize=(50, 22))  # (14, 20)
ax = plt.subplot(1, 1, 1)

    
AddTime = 2000*12  # spine up of 2000 years


                                                                                    
count = 1
ChosenStart = -100*12 # comparison between 2000 and 2100

# Regional coral carrying capacity
K_C_GBR = 438718.902336e10 # in cm^2
K_C_SEA = 1191704.68901e10 # in cm^2
K_C_CAR = 190374.697987e10  # in cm^2
K_C_List = [K_C_GBR, K_C_SEA, K_C_CAR]

# + or - 25% change
ParamsPlus= ["Gmax-Plus", "a-Plus", "b-Plus", "KC-Plus", "Ksmax-Plus", "MC-Plus", "alpha-Plus", "r-Plus", "beta-Plus", "GammaH-Plus"]
ParamsMinus= ["Gmax-Minus", "a-Minus", "b-Minus", "KC-Minus", "Ksmax-Minus", "MC-Minus", "alpha-Minus", "r-Minus", "beta-Minus", "GammaH-Minus"]
Params = ["$G_{max}$", "$a$", "$b$", "$K_C$", "$K_{smax}$", "$M_C$", r"$\alpha$", "$r$", r"$\beta$", r"$\Gamma_h$"]

NumP = linspace(1.5, 30, len(Params) + 1) 
col = [(1, 0.667, 0.667), (0.831, 0.416, 0.416), (0.502, 0.1, 0.1)]
ChosenStart = -100*12 # comparing betweeen year 2000 and 2100

cmx = mpl.cm
cmapC = mpl.cm.Set1
#cNorm = mpl.colors.Normalize(vmin=N_List.min(), vmax=N_List.max())
cNorm = mpl.colors.Normalize(vmin=0, vmax=len(Params))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmapC)
# Sensitivity of results with respect to speed of adaptation N, 
rawNum_GBR_orig = arange(0.01, 0.1, 0.00005)
rawNum_SEA_orig = arange(0.01, 0.1, 0.00005)
rawNum_CAR_orig = arange(0.01, 0.1, 0.00005) 

# These are not the index of estimated N reported in manuscript, it was from earlier code versions, the choice does not matter
GBR_N_index = 100 # index in rawNum_GBR  0.0145
SEA_N_index = 275 # index in rawNum_SEA  0.02375
CAR_N_index = 250 # index in rawNum_CAR  0.0225

# Only for sensitivity simulations to reduce the number of simulations (USED for supplementary figure in manuscript)
rawNum_GBR = concatenate((arange(min(rawNum_GBR_orig), rawNum_GBR_orig[GBR_N_index], 0.0005), array([0.75*rawNum_GBR_orig[GBR_N_index]]), array([rawNum_GBR_orig[GBR_N_index]]), array([1.25*rawNum_GBR_orig[GBR_N_index]]), arange(rawNum_GBR_orig[GBR_N_index]+0.0005, max(rawNum_GBR_orig), 0.0005)))
rawNum_SEA = concatenate((arange(min(rawNum_SEA_orig), rawNum_SEA_orig[SEA_N_index], 0.0005), array([0.75*rawNum_SEA_orig[SEA_N_index]]), array([rawNum_SEA_orig[SEA_N_index]]), array([1.25*rawNum_SEA_orig[SEA_N_index]]), arange(rawNum_SEA_orig[SEA_N_index]+0.0005, max(rawNum_SEA_orig), 0.0005)))
rawNum_CAR = concatenate((arange(min(rawNum_CAR_orig), rawNum_CAR_orig[CAR_N_index], 0.0005), array([0.75*rawNum_CAR_orig[CAR_N_index]]), array([rawNum_CAR_orig[CAR_N_index]]), array([1.25*rawNum_CAR_orig[CAR_N_index]]), arange(rawNum_CAR_orig[CAR_N_index]+0.0005, max(rawNum_CAR_orig), 0.0005))) 

"""
# Only for sensitivity simulations to reduce the number of simulations (Just for a DEMO here)
rawNum_GBR = concatenate((arange(min(rawNum_GBR_orig), rawNum_GBR_orig[GBR_N_index], 0.05), array([rawNum_GBR_orig[GBR_N_index]]), arange(rawNum_GBR_orig[GBR_N_index]+0.05, max(rawNum_GBR_orig), 0.05)))
rawNum_SEA = concatenate((arange(min(rawNum_SEA_orig), rawNum_SEA_orig[SEA_N_index], 0.05), array([rawNum_SEA_orig[SEA_N_index]]), arange(rawNum_SEA_orig[SEA_N_index]+0.05, max(rawNum_SEA_orig), 0.05)))
rawNum_CAR = concatenate((arange(min(rawNum_CAR_orig), rawNum_CAR_orig[CAR_N_index], 0.05), array([rawNum_CAR_orig[CAR_N_index]]), arange(rawNum_CAR_orig[CAR_N_index]+0.05, max(rawNum_CAR_orig), 0.05))) 
"""
# Parameters to study sensitivity to
ParamsPlus= ["Gmax-Plus", "a-Plus", "b-Plus", "KC-Plus", "Ksmax-Plus", "MC-Plus", "alpha-Plus", "r-Plus", "beta-Plus", "GammaH-Plus"]
ParamsMinus= ["Gmax-Minus", "a-Minus", "b-Minus", "KC-Minus", "Ksmax-Minus", "MC-Minus", "alpha-Minus", "r-Minus", "beta-Minus", "GammaH-Minus"]

#  Speed of acclimation specific for each region for main results, to indicate on the plot
GBR_N_index = 908 # index in rawNum_GBR  0.0554
SEA_N_index = 330 # index in rawNum_SEA  0.0265
CAR_N_index = 275 # index in rawNum_CAR  0.02375
reg_N_index = [GBR_N_index, SEA_N_index, CAR_N_index]

filename = "Sensitivity-N/Monthly/" #Specific folder for the sensitivity simulations "Sensitivity-RCP-monthly-2.py"

PlusOrMinus = "Plus"

if PlusOrMinus == "Plus":
    Params_which = ParamsPlus
else:
    Params_which = ParamsMinus

import pdb
from scipy import interpolate   
for z in xrange(len(Locations)):
    for v in xrange(len(RCP)):            
        rcp = RCP[v] 
        sub = plt.subplot(3, 3, count)
        for ind in xrange(len(Params_which)):
            CORval = array([])
            N_full = array([])
            
            if PlusOrMinus == "Plus":
                change = 1.25 # 0.75 for ParamsMinus and 1.25 for ParamsPlus
                fileLab = ParamsPlus[ind]
            else:
                change = 0.75
                fileLab = ParamsMinus[ind]
                
            file2 = open(filename+"CORAL-"+rcp+"-"+Locations[z]+"-"+fileLab+".dat", "r")
            HOSTSet1 = load(file2, allow_pickle = True)
            file2.close()
        
            file3 = open(filename+"TRAIT-"+rcp+"-"+Locations[z]+"-"+fileLab+".dat", "r")
            TRAITSet1 = load(file3, allow_pickle = True)
            file3.close()
        
            file4 =  open(filename+"SYMB-"+rcp+"-"+Locations[z]+"-"+fileLab+".dat", "r")
            SYMBSet1 = load(file4, allow_pickle = True)
            file4.close()
            
            HOST = HOSTSet1
            SYMB = SYMBSet1
            TRAIT = TRAITSet1
            K_C = K_C_List[z]

            if z == 0:
                rawNum_orig = rawNum_GBR_orig
                rawNum = rawNum_GBR
                N_List_orig = (scale)*rawNum_orig
                N_List = (scale)*rawNum
                plt.plot(N_List_orig[reg_N_index[z]]*ones(2), linspace(0, 1, 2), color = Color_list[0], linewidth = 8)
            elif z == 1:
                rawNum_orig = rawNum_SEA_orig
                rawNum = rawNum_SEA
                N_List_orig = (scale)*rawNum_orig
                N_List = (scale)*rawNum
                plt.plot(N_List_orig[reg_N_index[z]]*ones(2), linspace(0, 1, 2), color = Color_list[0], linewidth = 8)
            elif z == 2:
                rawNum_orig = rawNum_CAR_orig
                rawNum = rawNum_CAR
                N_List_orig = (scale)*rawNum_orig
                N_List = (scale)*rawNum
                plt.plot(N_List_orig[reg_N_index[z]]*ones(2), linspace(0, 1, 2), color = Color_list[0], linewidth = 8)
            print rawNum_orig[reg_N_index[z]], N_List_orig[reg_N_index[z]]

            ResCoral = zeros((len(Params), len(rawNum)))
            ResSymb = zeros((len(Params), len(rawNum)))
            ResTrait = zeros((len(Params), len(rawNum)))
        
            LocResCoral = zeros(len(rawNum))
            LocResSymb = zeros(len(rawNum))
            LocResTrait = zeros(len(rawNum))
            for i in xrange(len(rawNum)): 
                if Params[ind] == "$K_C$":
                    LocResCoral[i] = mean(HOST[i][ChosenStart:])/(change*K_C)
                    LocResSymb[i] = mean(TRAIT[i][ChosenStart:])
                    LocResTrait[i] = mean(SYMB[i][ChosenStart:])  
                else:
                    LocResCoral[i] = mean(HOST[i][ChosenStart:])/K_C
                    LocResSymb[i] = mean(TRAIT[i][ChosenStart:])
                    LocResTrait[i] = mean(SYMB[i][ChosenStart:])                    
            ResCoral[ind, :] = LocResCoral
            ResSymb[ind, :] = LocResSymb
            ResTrait[ind, :] = LocResTrait
            
            CORval = concatenate((CORval, ResCoral[ind, :]*((ind+3+1)/(len(Params)+3))))
            N_full = concatenate((N_full, N_List))
            
            if Params[ind] != r"$\Gamma_h$":
                colorval0 = scalarMap.to_rgba(ind+1)
            else:
                colorval0 = "black" 
        
            # cubic spline
            f = interpolate.interp1d(array(N_full), array(CORval), kind="cubic")
            N_new = linspace(min(N_full), max(N_full), 200)
            CORval_new = f(N_new)
            
            if count == 1:
                sub.plot(N_new, CORval_new, linewidth = 6, color = colorval0, label=Params[ind])#+ " (scaling = %.2f/$K_C$)"%((ind+3+1)/(len(Params)+3), ))
            else:
                sub.plot(N_new, CORval_new, linewidth = 6, color = colorval0)  
                
        if count not in (7, 8, 9):
            plt.xticks(array([1, 2, 3, 4, 5, 6, 7])*1e-13, [" " for i in [1, 2, 3, 4, 5, 6, 7]])
        else: 
            plt.xticks(array([1, 2, 3, 4, 5, 6, 7])*1e-13, ["%d"%N for N in [1, 2, 3, 4, 5, 6, 7]], fontsize=fsize2)
            plt.xlabel(r'Speed of acclimation ($\times$ $10^{-13}$)', fontsize = fsize2)
                
        if count in (1, 4, 7):
            if z == 0:
                plt.text(min(N_List) - 0.25e-12, 0.52, Locations_title[z], fontsize = fsize, fontproperties = font,  rotation = "vertical")
            if z == 1:
                plt.text(min(N_List) - 0.25e-12, 0.42, Locations_title[z], fontsize = fsize, fontproperties = font,  rotation = "vertical")
            if z == 2:
                plt.text(min(N_List) - 0.25e-12, 0.35, Locations_title[z], fontsize = fsize, fontproperties = font,  rotation = "vertical")    
            sub.text(min(N_List) - 0.1e-12, 0.5, Fig_lab[count-1], fontsize = fsize, fontproperties = font)
            hostTicks = array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
            plt.ylabel("Scaled mean\n coral abundance \n", horizontalalignment="center", fontsize = fsize2)
        else:
            plt.text(min(N_List) - 0.035e-12, 0.5, Fig_lab[count-1], fontsize = fsize, fontproperties = font)
            plt.yticks(list(hostTicks), [" "]+list([" " for k in xrange(1, len(hostTicks))]), fontsize = fsize2)
        if count in (1, 2, 3):
            plt.title(RCP_title[v], fontsize = fsize, fontproperties = font)
        xmax = max(N_List)
        if count in (1, 4, 7):
            sub.tick_params(axis = "y", direction = "in", pad = 2, labelsize = fsize2)
            sub.set_yticks(hostTicks)
            sub.set_yticklabels(["$0$"]+list(["%.1f"%hostTicks[k] for k in xrange(1, len(hostTicks))]))
            part = sub.twinx()
            part.tick_params(axis = "y", direction = "in", pad = 2, labelsize = fsize2)
            part.set_yticks(hostTicks)
            part.set_yticklabels(["" for h in hostTicks])
            part.set_ylim((0, 0.5))
        else:
            sub.tick_params(axis = "y", direction = "in", pad = 2, labelsize = fsize2)
            sub.set_yticks(hostTicks)
            sub.set_yticklabels(["" for h in hostTicks])
            part = sub.twinx()
            part.tick_params(axis = "y", direction = "in", pad = 2, labelsize = fsize2)
            part.set_yticks(hostTicks)
            part.set_yticklabels(["" for h in hostTicks])
            part.set_ylim((0, 0.5))
        sub.set_xlim((1e-13, 7e-13))
        sub.set_ylim((0, 0.5))
        if count == 1:    
            #sub.legend(loc=(-0.75, 1.25), fontsize = fsize2, ncol = 5)
            sub.legend(loc=(0.6, 1.25), fontsize = fsize, ncol = 5) 
      
        count+=1

plt.subplots_adjust(top=0.965,
bottom=0.095,
left=0.080,
#right=0.81,
right = 0.75,
hspace=0.175,
wspace=0.175)

if PlusOrMinus == "Plus":
    fig.savefig("Figures/EPS/Fig5.eps", dpi= 600, bbox_inches = 'tight') 
    fig.savefig("Figures/PDF/Fig5.pdf", dpi= 600, bbox_inches = 'tight')
else:
    fig.savefig("Figures/EPS/Fig6.eps", dpi= 600, bbox_inches = 'tight') 
    fig.savefig("Figures/PDF/Fig6.pdf", dpi= 600, bbox_inches = 'tight')                                   

                                                                                                                                                                           
#plt.show()