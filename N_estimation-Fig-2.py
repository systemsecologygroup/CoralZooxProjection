# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, arange, isnan
from scipy.integrate import ode
from matplotlib import pyplot as plt
import matplotlib
import numpy as np

""" 
This code was used to estimate the speed of acclimation N, Figure 2 in Manuscript

All results must first be generated for all the N values in Helper.py and then stored in the folder bellow e.g. Results_final 
(The files are too heavy to store on github)
"""

filename = "Results-Monthly-fullRange-N/" # Contain simulations for all the values of N in the full ranges in Helper.py 


#### Model with a dynamics on the symbiont biomass, Time in month #####
G_C = 10/12  # this is G_max in the model
a = 1.0768/12 # Symbiont specific growth rate - linear growth rate
b = 0.0633 # Exponential Growth constant of the symbiont


alpha = 1e-3

M_C = (10/12)*1e-3 # coral mortality 
# Regional coral carrying capacity
K_C_GBR = 438718.902336e10 # in cm^2
K_C_SEA = 1191704.68901e10 # in cm^2
K_C_CAR = 190374.697987e10  # in cm^2

K_C_List = [K_C_GBR, K_C_SEA, K_C_CAR]

Ksmax = 3e6 # healthy measure of carying capacity of symbiont per host biomass Gamma

beta = (12)*1e2

gammaH = 0.25e6 # free param 
fsize = 18

r = (12)*1e3

color_list = array([(0.647*col, 0.333*col, 0.075*col) for col in xrange(2,5+2)])/5

scale = 1e-11

rawNum_GBR = np.arange(0.01, 0.1, 0.00005)
rawNum_SEA = np.arange(0.01, 0.1, 0.00005)
rawNum_CAR = np.arange(0.01, 0.1, 0.00005) 

RCP_list = ["RCP26"]#, "RCP45", "RCP85"]

Locations = array(["GBR", "SEA", "CAR"])


Locations = ["GBR", "SEA", "CAR"]
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]
RCP = ["RCP26", "RCP45", "RCP85"] # for filename
RCP_title =  ["RCP2.6", "RCP4.5", "RCP8.5"]
Fig_lab = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
#Fig_lab = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]

# for a nice plot
fsize = 16*3 #23
fsize2 = 14*3 #14 #18

import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

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
      

#Color_list = [(0.306, 0.431+0.1, 0.545+0.3), (0.839+0.025, 0.60+0.05, 0.20), (0.839+0.16, 0.363, 0.35)]
Color_list = [(0.651, 0.325, 0.529), (0.839+0.025, 0.60+0.05, 0.20), (0.839+0.16, 0.363, 0.35)]

# Determine slope and intercept from Fitting-coral-cover.r 

# Data for GBR from source1: De'ath et al. 2012
source1Col = "#006024"#(0.173, 0.375, 0.187) 
GBR_Year_source1 = array([1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                            1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012])
GBR_Cover_source1 = array([27.6392, 22.7458, 26.3999, 20.5105, 19.6836, 21.1800, 23.0085, 23.4262, 22.1840, 23.1825, 25.7579, 23.1051, 
                            20.8673, 25.4344, 24.2754, 23.0333, 20.6293,17.9764,17.4814, 21.7995, 16.3249,20.4771,15.5005, 18.5738, 
                            18.3276, 16.0067, 10.3663, 9.70535])    

# Data for GBR from source2: Bruno 2007, Regional Decline of Coral Cover in the Indo-Pacific Timing, Extent, and Subregional Comparisons

source2Col = (0.353, 0.675, 0.337) #(0.063, 0.302, 0.)
GBR_Year_source2 = array([1971, 1972, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 
                    1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007])
GBR_Cover_source2 = array([4.2999999999999998, 34.975000000000001, 49.375, 24.906481481481478, 24.441163003663007, 
                        24.126298701298698, 16.893939393939394, 22.58467261904762,17.311111111111114, 21.028985714285714, 
                        23.527606393606391, 25.455069209922868, 28.285919266492023, 31.795480366220563, 32.272448584582911,
                        32.982510303536621, 33.638900123685843, 29.743189587198515, 30.616005555555553, 27.787401217829789,
                        31.235181694495974, 30.107761904761908, 26.650676470588238, 30.476910112359551,14.855399999999999])     
GBR_slope = 0.1848
GBR_intercept = -341.6716
GBR_confint = [-0.1939096, 0.5634574]  #95%                  
                                              
# Data for SEA from source2: Bruno 2007, Regional Decline of Coral Cover in the Indo-Pacific Timing, Extent, and Subregional Comparisons
SEA_Year_source2 = array([1971, 1974, 1977, 1978, 1979, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 
                    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006])                
SEA_Cover_source2 = array([65.0, 39.683333333333337, 15.176470588235293, 10.75, 40.615229885057467, 38.426666666666669, 47.682051282051283,
                        29.899999999999999, 34.222222222222221, 27.421666666666663, 38.978494623655912, 42.787939221272559, 35.265217391304347,
                        30.222222222222221, 49.857142857142854, 38.984946236559139, 30.648148148148149, 35.902311594202892, 39.000386904761903, 
                        38.747898397435897, 27.335411368015414, 34.826259328358212, 27.322793103448276, 29.019946808510639, 30.414325842696631,
                        33.062789351851855, 31.652665770609318, 34.652255639097746, 31.402439024390244, 27.543560606060606, 33.609910337552741])
SEA_slope = -0.2354
SEA_intercept = 503.0231
SEA_confint = [-0.6052167, 0.1345118]  

# Data for CAR from source2: Bruno 2007, Regional Decline of Coral Cover in the Indo-Pacific Timing, Extent, and Subregional Comparisons          
CAR_Year_source2 = array([1976, 1977, 1978, 1980, 1981, 1982, 1983, 1984, 1986, 1987, 1988, 1989,1990, 1992,
                    1993, 1994, 1995,1996,1997,1998,1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006])
CAR_Cover_source2 = array([25.17, 47.1, 46.7735, 28.0875, 16.1, 22.4429375,  17.2,  12.26666667,
         6.78571429, 10.296, 22.4, 2.55, 27.85,   6.3375,  22.32375, 24.10105556,
        26.1696, 17.36739394, 25.76126312, 32.40814444, 28.30322, 15.62805952, 20.24744444, 35.37140909,
        26.1427, 24.74643878,21.81088542, 14.9610625])
        
CAR_slope = -0.1304
CAR_intercept = 282.1811
CAR_confint = [-0.5865406, 0.3256516] 


# Save data to make available publicly

CovGBR = array((GBR_Year_source2, GBR_Cover_source2))
CovSEA = array((SEA_Year_source2, SEA_Cover_source2))
CovCAR = array((CAR_Year_source2, CAR_Cover_source2))

FileGBR = open("Percent_Coral_Cover_GBR_between_1971_2007.dat", "wr")
FileSEA = open("Percent_Coral_Cover_SEA_between_1971_2006.dat", "wr")
FileCAR = open("Percent_Coral_Cover_CAR_between_1976_2006.dat", "wr")

CovGBR.dump(FileGBR)
CovSEA.dump(FileSEA)
CovCAR.dump(FileCAR)

FileGBR.flush()
FileSEA.flush()
FileCAR.flush()

count = 0

row = 1

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


IndReg = [(700, 900), (200, 300), (200, 300)]
IndReg2 = [(900, 1000), (300, 400), (300, 400)]
IndReg3 = [(1000, 1200), (400, 600), (400, 600)]
IndReg4 = [(1200, len(rawNum_GBR)), (600, len(rawNum_SEA)), (600, len(rawNum_CAR))]

LocReg = [15, 10, 10]
LocReg2 = [15, 25, 15]
LocReg3 = [90, 90, 90]#[20, 30, 30]
LocReg4 = [120, 140, 140]

IndReg_List = [IndReg, IndReg2, IndReg3, IndReg4]
LocReg_List = [LocReg, LocReg2, LocReg3, LocReg4]

#  Speed of adapations specific for each region # index in N_list
#in case of change then redo the runs for SuppFig 10 and 11: sensitivity wrt to the chosen N
GBR_N_index = 908 # index in rawNum_GBR  0.0554 
SEA_N_index = 330 # index in rawNum_SEA  0.0265
CAR_N_index = 275 # index in rawNum_CAR  0.02375
reg_N_index = [GBR_N_index, SEA_N_index, CAR_N_index]
# Slope and intercept list for model result (by region), first determine the indexes above and then copy and past the results in Fitting-coral-cover.r to get the intercept and slope bellow
#intercept_model = [125.23889, 496.633902, 259.10463] # SEA_N_index = 340
#slope_model = [-0.04985, -0.231949, -0.11867] # SEA_N_index = 340

intercept_model = [125.23889, 492.368581, 259.10463]
slope_model = [-0.04985, -0.230158, -0.11867]
tcks = (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)

# Parameter resp. for the different region GBR, SEA , CAR
T_opt = array([26.76, 28.46, 27.56])

AddTime = 2000*12

fig = plt.figure(figsize=(60, 30))
markS = 20 # green dots size
import pdb
for z in xrange(len(Locations)):
    part0 = plt.subplot(row, 3, count+1)
    
    plt.title(Locations_title[z], fontsize = fsize, fontproperties = font)#fsize)
    if z == 0:
        plt.text(1964, 110, Fig_lab[z], fontsize = fsize, fontproperties = font)
    else:
        plt.text(1967, 110, Fig_lab[z], fontsize = fsize, fontproperties = font)
    part1 = part0.twinx()
    
    for rcp_index in xrange(1): # we only need to plot for one scenario because we plot for a time-range that doesn't include future projection
        rcp = RCP[rcp_index]
        v = rcp_index  
                
        file0 = open("Monthly-SST-scenarios/Months-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        time0 = load(file0,allow_pickle = True)
        file0.close()
        
        time = concatenate((arange(min(time0)-(AddTime+12)/12, min(time0), 1/12), time0))   
        
        file1sst = open("Monthly-SST-scenarios/SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        SST = load(file1sst,allow_pickle = True)
        file1sst.close()
        
        file2 =open(filename+rcp+"/CORAL-"+rcp+"-"+Locations[z]+".dat", "r")
        HOSTSet1 = load(file2,allow_pickle = True)
        file2.close()
        
        file3 = open(filename+rcp+"/TRAIT-"+rcp+"-"+Locations[z]+".dat", "r")
        TRAITSet1 = load(file3,allow_pickle = True)
        file3.close()
        
        file4 =  open(filename+rcp+"/SYMB-"+rcp+"-"+Locations[z]+".dat", "r")
        SYMBSet1 = load(file4,allow_pickle = True)
        file4.close()
        
        HOST = HOSTSet1 
        SYMB = SYMBSet1
        TRAIT = TRAITSet1        

        K_C_Reg = K_C_List[z]
        
        if z == 0:
            N_List = (scale)*rawNum_GBR
        elif z == 1:
            N_List = (scale)*rawNum_SEA
        else:
            N_List = (scale)*rawNum_CAR  
            part1.set_ylabel("Simulated relative \n coral abundance ($\%$)", color = Color_list[v], fontproperties = font, fontsize = fsize+2, labelpad = 2)

        for i in xrange(0, len(N_List)):
            for indr in xrange(len(IndReg_List)):
                IndREG = IndReg_List[indr]
                LocREG = LocReg_List[indr]
                if i >= IndREG[z][0] and i <=IndREG[z][1] and i%LocREG[z] == 0:
                    #print i
                    Host = HOST[i]
                    Trait = TRAIT[i]
                    Symb = SYMB[i]                                  
                    part1.plot(time, 100*Host/K_C_Reg, linewidth = 0.75, color = Color_list[v], alpha = 0.35)    
                    

        # Plot for specific speed of adaptation N
        Host = HOST[reg_N_index[z]]
        print Locations[z], reg_N_index[z], time[(time >= 1970)*(time<=2010)][::12]
        print Locations[z], reg_N_index[z], 100*Host[(time >= 1970)*(time<=2010)][::12]/K_C_Reg
        if z==0:
            #fit_model = slope_model[z]*time[(time >= 1970)*(time<=2010)]+ intercept_model[z]
            part1.plot(time, 100*Host/K_C_Reg, linewidth = 15, color = Color_list[v], alpha = 0.75)
            part1.plot(time0, 100*SST/T_opt[z],linewidth = 3, alpha= 0.4, color=(0.7, 0.7, 0.7))
            #part1.plot(time[(time >= 1970)*(time<=2010)], fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75)
            #part0.plot(0*time[(time >= 1970)*(time<=2010)], 0*fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75, label="$y$ = %.2f$x$ + %.2f"%(slope_model[z], intercept_model[z]))
            #part0.legend(frameon = True, fontsize = 16)
            part1.set_yticks(tcks) 
            part1.set_yticklabels(["" for i in range(len(tcks))], fontsize = fsize)  
            part0.plot(time, -(0.5)*100*Host/K_C_Reg, linewidth = 15, color = Color_list[v], alpha = 0.75, label="Best simulation") # this is only for showing label, it does not appear in figure
            Host_i = HOST[i]
            Trait_i = TRAIT[i]
            Symb_i = SYMB[i] 
            part0.plot(time, -(0.5)*Host_i/K_C_Reg, linewidth = 0.75, color = Color_list[v], alpha = 0.35 , label="Simulations")## this is only for showing label,it does not appear in figure
            part0.plot(time0, -(0.5)*SST/T_opt[z],linewidth = 3, alpha= 0.4, color=(0.7, 0.7, 0.7), label = "% SST$/T^{opt}$") # this is only for showing label,it does not appear in figure
        elif z == 1:
            #fit_model = slope_model[z]*time[(time >= 1970)*(time<=2010)]+ intercept_model[z]
            part1.plot(time, 100*Host/K_C_Reg, linewidth = 15, color = Color_list[v], alpha = 0.75)
            part1.plot(time0, 100*SST/T_opt[z], linewidth = 3, alpha= 0.4, color=(0.7, 0.7, 0.7))
            #part1.plot(time[(time >= 1970)*(time<=2010)], fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75)
            #part0.plot(0*time[(time >= 1970)*(time<=2010)], 0*fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75, label="$y$ = %.2f$x$ + %.2f"%(slope_model[z], intercept_model[z]))
            part1.set_yticks(tcks) 
            part1.set_yticklabels(["" for i in range(len(tcks))], fontsize = fsize)    
        elif z==2:
            fit_model = slope_model[z]*time[(time >= 1970)*(time<=2010)]+ intercept_model[z]
            part1.plot(time, 100*Host/K_C_Reg, linewidth = 15, color = Color_list[v], alpha = 0.75)
            part1.plot(time0, 100*SST/T_opt[z],linewidth = 2, alpha= 0.4, color=(0.7, 0.7, 0.7))
            #part1.plot(time[(time >= 1970)*(time<=2010)], fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75)
            #part0.plot(0*time[(time >= 1970)*(time<=2010)], 0*fit_model, linewidth = 6, color = Color_list[v], alpha = 0.75, label="$y$ = %.3f$x$ + %.2f"%(slope_model[z], intercept_model[z]))
            #part0.legend(frameon = True, fontsize = fsize2, facecolor = "white")
            
            part1.set_yticks(tcks) 
            part1.set_yticklabels(tcks, color = Color_list[v], fontsize = fsize)                                    
    
    # Plot data or fit
    if z == 0: # GBR
        part0.set_yticks(tcks) 
        part0.set_yticklabels(tcks, color = source2Col, fontsize = fsize)
        part0.plot(GBR_Year_source2, GBR_Cover_source2, "o", markersize = markS, alpha = 0.75, color = source2Col, label="Obs. (Bruno & Selig, 2007)")
        part0.plot(GBR_Year_source1, GBR_Cover_source1, "o", markersize = markS, alpha = 0.75, color = source1Col, label="Obs. (De'ath $et$ $al.$, 2012)")
        #fit_GBR = GBR_slope*GBR_Year_source2 + GBR_intercept
        #part0.plot(GBR_Year_source2, fit_GBR, linewidth = 6, color = source2Col, alpha = 0.75, label="$y$ = %.2f$x$ - %.2f"%(GBR_slope, abs(GBR_intercept)))
        #part0.legend(frameon = True, fontsize = fsize2, loc = (0.02,0.5), facecolor = "white", framealpha = 1)
        part0.legend(frameon = True, fontsize = fsize2+2, loc = (0.2,-0.3), ncol = 5, facecolor = "white", framealpha = 1)
        part0.set_ylabel("Observed relative \n coral abundance ($\%$)", color = source2Col, fontproperties = font, fontsize = fsize+2, labelpad = 2)
    elif z == 1: # SEA
        part0.set_yticks(tcks) 
        part0.set_yticklabels(["" for i in range(len(tcks))], color = source2Col, fontsize = fsize)
        #fit_SEA = SEA_slope*SEA_Year_source2 + SEA_intercept
        #part0.plot(SEA_Year_source2, fit_SEA, linewidth = 6, color = source2Col, alpha = 0.75, label="$y$ = %.2f$x$ + %.2f"%(SEA_slope, SEA_intercept))
        part0.plot(SEA_Year_source2, SEA_Cover_source2, "o", markersize = markS, alpha = 0.75, color = source2Col)#, label="Data ($\%$ Cover)")
        #part0.legend(frameon = False, fontsize = 16)

    elif z == 2: # CAR
        part0.set_yticks(tcks) 
        part0.set_yticklabels(["" for i in range(len(tcks))], color = source2Col, fontsize = fsize)
        part0.plot(CAR_Year_source2, CAR_Cover_source2, "o", markersize = markS, alpha = 0.75, color = source2Col)#, label="Data ($\%$ Cover)")

        #fit_CAR = CAR_slope*CAR_Year_source2 + CAR_intercept
        #part0.plot(CAR_Year_source2, fit_CAR, linewidth = 6, color = source2Col, alpha = 0.75, label="$y$ = %.3f$x$ + %.2f"%(CAR_slope, CAR_intercept))

    part0.set_xlim((1970, 2010))
    part0.set_xticks(arange(1970, 2015, 5))
    part0.set_xticklabels(["%d"%k for k in arange(1970, 2015, 5)], rotation = 45, fontsize = fsize)#fontsize = fsize)
    part0.set_ylim((0,110)) 
    part1.set_ylim((0,110))
    part0.set_xlabel("Years", fontsize = fsize+2)
    
    if count == 0: 
        plt.arrow(2008, 48, 0, 10, head_width = 1, head_length = 2, color = "black") 
    elif count == 1:
        plt.arrow(2008, 56, 0, 10, head_width = 1, head_length = 2, color = "black") 
    elif count == 2:
        plt.arrow(2008, 56, 0, 10, head_width = 1, head_length = 2, color = "black") 
    count +=1

# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

# Adjust margins
plt.subplots_adjust(bottom = 0.47, right = 0.91, left = 0.10, top = 0.91, wspace = 0.14, hspace = 0.20)

plt.savefig("Figures/EPS/Fig2.eps", dpi= 600, bbox_inches = 'tight') 
plt.savefig("Figures/PDF/Fig2.pdf", dpi= 600, bbox_inches = 'tight')

#plt.show()

