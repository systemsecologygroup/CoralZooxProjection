from __future__ import division
from scipy import array, load, arange, linspace, ones, cumsum, concatenate, mean
from scipy.stats import norm

from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib
plt.switch_backend('agg')


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

fsize = (16+8) #18 #22
fsize2 = (16+8) #18

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
 
fig =plt.figure(figsize=(30, 13+2))  # (14, 20)
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
for z in range(len(Locations)):
    sub1 = plt.subplot(4, 3, count+1)
    part1 = sub1.twinx()
    plt.title(Locations_title[z], fontsize = fsize+2, fontproperties = font)
    T0 = T0_list[z]
    rho = rho_list[z]
    skew = skew_list[z]
    TempCORListCenter = (TempList - T0)/rho
    NormCor = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)
    Dist = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)/max(NormCor)
    maxD = max(Dist)
    T_opt = TempList[int(list(Dist).index(maxD))]
    for v in range(len(RCP)):
        rcp = RCP[v]
        T0 = T_opt # T_opt is used as a threshold for bleaching events
        file0 = open("Monthly-SST-scenarios/Months-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        time0 = load(file0, allow_pickle = True)
        file0.close()
        
        file1sst = open("Monthly-SST-scenarios/SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        SST = load(file1sst, allow_pickle = True)
        file1sst.close()
        if z == 0:
            CummulativeBl1 = cumsum((SST[list(time0).index(startTime):] >= T0+2)*(SST[list(time0).index(startTime):] < T0+4))    # Cummulaticve number of bleaching since 2016
            CummulativeBl2 = cumsum(SST[list(time0).index(startTime):] >= T0+4)    # Cummulaticve number of bleaching since 2016
            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl1 != 0], CummulativeBl1[CummulativeBl1 != 0], linewidth = 3, color = Color_list[v], label=RCP_title[v])#+r"($2 \leq \Delta T_{model}\,<\,4$)")
            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl2 != 0], CummulativeBl2[CummulativeBl2 != 0], "--", linewidth = 3, color = Color_list[v])#, label=RCP_title[v]+r"($\Delta T_{model} \geq 4$)")
            sub1.plot(time0[list(time0).index(startTime):], CummulativeBl1 + CummulativeBl2, linewidth = 4, color = Color_list[v], label=RCP_title[v])#+r"($2 \leq \Delta T_{model}\,<\,4$)")
        elif z == 1:
            CummulativeBl = cumsum(SST[list(time0).index(startTime):] >= T0+1)    # Cummulative number of bleaching since 2016
            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl != 0], CummulativeBl[CummulativeBl != 0], linewidth = 3, color = Color_list[v])# label=RCP_title[v]+r"($\Delta T_{model} \geq 1$)")
            sub1.plot(time0[list(time0).index(startTime):], CummulativeBl, linewidth = 4, color = Color_list[v])# label=RCP_title[v]+r"($\Delta T_{model} \geq 1$)")
        else:
            CummulativeBl1 = cumsum((SST[list(time0).index(startTime):] >= T0+1)*(SST[list(time0).index(startTime):] < T0+3))    # Cummulaticve number of bleaching since 2016
            CummulativeBl2 = cumsum((SST[list(time0).index(startTime):] >= T0+3)*(SST[list(time0).index(startTime):] < T0+6))    # Cummulaticve number of bleaching since 2016
            CummulativeBl3 = cumsum(SST[list(time0).index(startTime):] >= T0+6)    # Cummulaticve number of bleaching since 2016

            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl1 != 0], CummulativeBl1[CummulativeBl1 != 0], linewidth = 3, color = Color_list[v])#, label=RCP_title[v]+r"($1 \leq \Delta T_{model}\,<\,3$)")
            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl2 != 0], CummulativeBl2[CummulativeBl2 != 0], "--",linewidth = 3, color = Color_list[v])#, label=RCP_title[v])#+ r"($3 \leq \Delta T_{model}\,<\,6$)")
            #sub1.plot(time0[list(time0).index(startTime):][CummulativeBl3 != 0], CummulativeBl3[CummulativeBl3 != 0], ".", markersize = 3 ,color = Color_list[v])#, label=RCP_title[v])#+r"($\Delta T_{model} \geq 6$)")
            sub1.plot(time0[list(time0).index(startTime):], CummulativeBl1+CummulativeBl2+CummulativeBl3, linewidth = 4, color = Color_list[v])#, label=RCP_title[v]+r"($1 \leq \Delta T_{model}\,<\,3$)")

        if count == 0:
            sub1.set_ylabel("Bleaching \n (nr)", fontsize = fsize)
            sub1.text(startTime - 16, 800, Fig_lab[count], fontproperties=font, fontsize = fsize)
            sub1.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub1.set_yticks((0, 200, 400, 600, 800, 1000))
            part1.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part1.set_yticks((0, 200, 400, 600, 800, 1000))
            part1.set_yticklabels(("", "", "", "", "", ""))
        elif count == 2:
            sub1.text(startTime - 5, 800, Fig_lab[count], fontproperties=font,  fontsize = fsize)
            sub1.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub1.set_yticks((0, 200, 400, 600, 800, 1000))
            sub1.set_yticklabels(("", "", "", "", "", ""))
            part1.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part1.set_yticks((0, 200, 400, 600, 800, 1000))
            part1.set_yticklabels(("", "", "", "", "", ""))
        else:
            sub1.text(startTime - 5, 800 , Fig_lab[count], fontproperties=font, fontsize = fsize)
            sub1.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub1.set_yticks((0, 200, 400, 600, 800, 1000))
            sub1.set_yticklabels(("", "", "", "", "", ""))
            part1.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part1.set_yticks((0, 200, 400, 600, 800, 1000))
            part1.set_yticklabels(("", "", "", "", "", ""))
    sub1.set_ylim((0, 800))
    part1.set_ylim((0, 800))
    plt.xticks([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100], ["", "", "", "", "", "", "", "", "", ""], fontsize = fsize)
    plt.xlim((startTime, max(time0)))
    if z == 0:
        sub1.legend(loc="upper left", fontsize = fsize2, frameon = True) 
    count +=1

count = 3

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


for z in range(len(Locations)):
    for v in range(len(RCP)):
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
        
        sub1 = plt.subplot(4, 3, 1+count)
        #TraitTicks = array([-0.5, 0, 101-100, 102-100, 102.5-100])
        TraitTicks = array([-0.5, 0, 0.5, 1, 1.5, 2, 2.5])
        plt.xticks([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100], ["", "", "", "", "", "", "", "", "", ""])
        sub1.set_xlim((startTime, max(time0))) 
        part1 = sub1.twinx()
        if count == 3:
            sub1.set_ylabel("Coral trait \n (% change)", fontsize = fsize)
            sub1.text(startTime - 16, 102.5-100, Fig_lab[count], fontproperties=font, fontsize = fsize) # forcing
        else:
            sub1.text(startTime - 5, 102.5-100, Fig_lab[count], fontproperties=font,  fontsize = fsize)
            
        sub4 = plt.subplot(4, 3, 4+count)
          
        plt.xticks([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100], ["", "", "", "", "", "", "", "", "", ""])
        part4 = sub4.twinx()
        #hostTicks = array([0, 20, 40, 60, 80, 100, 120]) - 100    
        hostTicks = array([0-100, 20-100, 40-100, 60-100, 80-100, 0, 120-100])  
        if count == 3:
            sub4.set_ylabel("Coral abundance \n (% change)", fontsize = fsize)
            sub4.text(startTime - 16, 120-100, Fig_lab[3+count],fontproperties=font, fontsize = fsize) # forcing

        else: 
            sub4.text(startTime - 5, 120-100, Fig_lab[3+count],fontproperties=font, fontsize = fsize) # forcing
    
        sub7 = plt.subplot(4, 3, 7+count)
        plt.xticks([2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100], rotation = 45, fontsize = fsize)

        Symbticks = array([0, 30, 60, 90, 120, 150, 180, 200]) - 100
        part7 = sub7.twinx()
        if count==3:
            sub7.set_ylabel("Symbiont \n abundance \n (% change)", fontsize = fsize)
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
        compare0 = (time<=2005)*(time>=1986)
        
        PastHost = sum(Host[compare0])/sum(compare0)
        PastTrait = sum(Trait[compare0])/sum(compare0)
        PastSymb = sum(Symb[compare0])/sum(compare0)
        
        sub1.plot(time, 100*Trait/PastTrait - 100, linewidth = 4, color = Color_list[v]) 
        sub4.plot(time, 100*Host/PastHost - 100, linewidth = 4, color = Color_list[v])
        sub7.plot(time, 100*Symb/PastSymb - 100, linewidth = 4, color = Color_list[v], label = RCP_title[v])
        
        compare1 = (time<=2100)*(time>=2081)     
        FutureHost = sum(Host[compare1])/sum(compare1)
        MeanlossCoralBiomass.append((FutureHost-PastHost)/PastHost)
        PerChange = mean(array(MeanlossCoralBiomass))*100
        print(rcp, Locations[z], PerChange, "N = ", rawNum_CAR[N_index_true]*scale)
        
        # saving time in .dat
        time_file = open("Results/Time-"+rcp+"-"+Locations[z]+"-%d.dat"%reg_N_index_true[z], "wr")
        time.dump(time_file)
        time_file.flush()
        time_file.close()
        
        # saving to csv
        """
        np.savetxt("Results-csv/Time-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], time, delimiter = ",")
        np.savetxt("Results-csv/CORAL-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Host, delimiter = ",")
        np.savetxt("Results-csv/TRAIT-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Trait, delimiter = ",")
        np.savetxt("Results-csv/SYMB-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Symb, delimiter = ",")
        """
        
        if count in (3, ):
            sub1.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub1.set_yticks(list(TraitTicks)) # GBR and SEA
            #sub1.set_yticklabels(["-0.5"]+list(["%d"%TraitTicks[s] for s in range(1, len(TraitTicks)-1)])+["2.5"])
            sub1.set_yticklabels(list(["%.1f"%TraitTicks[s] for s in range(len(TraitTicks))]))
            part1.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part1.set_yticks(TraitTicks)
            part1.set_yticklabels([" "%d for d in TraitTicks])
            
            sub4.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub4.set_yticks(list(hostTicks))
            sub4.set_yticklabels(list(["%d"%hostTicks[k] for k in range(len(hostTicks))]))
            part4.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part4.set_yticks(hostTicks)
            part4.set_yticklabels([" "%d for d in hostTicks])
            
            sub7.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub7.set_yticks(list(Symbticks))
            sub7.set_yticklabels(list(["%d"%Symbticks[s] for s in range(len(Symbticks))]))
            part7.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part7.set_yticks(Symbticks)
            part7.set_yticklabels([" "%d for d in Symbticks])
        else:
            sub1.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub1.set_yticks(list(TraitTicks))
            sub1.set_yticklabels([" "]+list([" " for s in range(1, len(TraitTicks))]))
            part1.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part1.set_yticks(TraitTicks)
            part1.set_yticklabels([" "%d for d in TraitTicks])
            
            sub4.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub4.set_yticks(list(hostTicks))
            sub4.set_yticklabels([" "]+list([" " for k in range(1, len(hostTicks))]))
            part4.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part4.set_yticks(hostTicks)
            part4.set_yticklabels([" "%d for d in hostTicks])
            
            sub7.tick_params(axis = "y", direction = "in", labelsize = fsize2)
            sub7.set_yticks(list(Symbticks))
            sub7.set_yticklabels([" "]+list([" " for s in range(1, len(Symbticks))]))
            part7.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
            part7.set_yticks(Symbticks)
            part7.set_yticklabels([" "%d for d in Symbticks])           
        
        sub1.set_ylim((-0.5, 102.5-100))
        sub1.set_xlim((startTime, max(time0))) 
        part1.set_ylim((-0.5, 102.5-100))
        part1.set_xlim((startTime, max(time0))) 
        
        sub4.set_xlim((startTime, max(time0))) 
        sub4.set_ylim((-100, 120-100))
        part4.set_ylim((-100, 120-100))
        sub4.set_xlim((startTime, max(time0))) 
        
        sub7.set_ylim((0-100, 180-100))
        part7.set_ylim((0-100, 180-100)) 
        sub7.set_xlim((startTime, max(time0))) 
        sub7.set_xlabel("Years", fontsize=fsize)  
        
        
        ### Plot simulation without adaptation
        
        if rcp == "RCP26":
            print("done")
             # plot "No adaptation from 2010"
             
            
            file2 = open("Results-N-0/Alltime/CORAL-"+rcp+"-"+Locations[z]+".dat", "r")
            HOSTSet1b = load(file2, allow_pickle = True)
            file2.close()
            
            file3 = open("Results-N-0/Alltime/TRAIT-"+rcp+"-"+Locations[z]+".dat", "r")
            TRAITSet1b = load(file3, allow_pickle = True)
            file3.close()
            
            file4 =  open("Results-N-0/Alltime/SYMB-"+rcp+"-"+Locations[z]+".dat", "r")
            SYMBSet1b = load(file4, allow_pickle = True)
            file4.close()
            
            HOST2 = HOSTSet1b
            SYMB2 = SYMBSet1b
            TRAIT2 = TRAITSet1b
            
            
            Host2 = HOST2[0]
            Trait2 = TRAIT2[0]
            Symb2 = SYMB2[0]
                        
            # save to csv
            #np.savetxt("Results-N-0-2010/Time-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], time, delimiter = ",")
            #np.savetxt("Results-N-0-2010/CORAL-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Host2, delimiter = ",")
            #np.savetxt("Results-N-0-2010/TRAIT-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Trait2, delimiter = ",")
            #np.savetxt("Results-N-0-2010/SYMB-"+rcp+"-"+Locations[z]+"-%d.csv"%reg_N_index_true[z], Symb2, delimiter = ",")
            
            #pdb.set_trace()
            
            PastHost2 = sum(Host2[compare0])/sum(compare0)
            PastTrait2 = sum(Trait2[compare0])/sum(compare0)
            PastSymb2 = sum(Symb2[compare0])/sum(compare0)
            
            sub1.plot(time, 100*Trait2/PastTrait2 - 100, linewidth = 2, color = "black", label = "%s (without acclimation)"%RCP_title[v])
            sub4.plot(time, 100*Host2/PastHost2 - 100, linewidth = 2, color = "black")
            sub7.plot(time, 100*Symb2/PastSymb2 - 100, linewidth = 2, color = "black")
            if z == 0:
                sub1.legend(fontsize = 15, loc="upper left")
            
    count +=1 


# Plot with maximal window
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()

plt.subplots_adjust(top=0.965,
bottom=0.095,
left=0.080,
#right=0.81,
right = 0.75,
hspace=0.175,
wspace=0.175)

plt.savefig("Figures/EPS/Fig3.eps", dpi= 600, bbox_inches = 'tight') 
plt.savefig("Figures/PDF/Fig3.pdf", dpi= 600, bbox_inches = 'tight')   
