from __future__ import division
from scipy import array, load, arange, linspace, ones
from scipy.stats import norm

from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties

"""
This code is used to generate Figure 4 in Manuscript

"""


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


Location = ["GBR", "SEA", "CAR"] 
LonCoordinates = [(145, 165), (100, 137), (360-80, 360-65)]
LatCoordinates = [(-28, -10), (-10, 13), (10, 20)]

######## Generating yearly temperature forcing by taking to year from the date of collection of WOD13 data #######             
### Locations in columns and RCPs in rows

Locations = array(["GBR", "SEA", "CAR"])
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]
RCP = ["RCP26", "RCP45", "RCP85"] # for filename
RCP_title =  ["RCP 2.6", "RCP 4.5", "RCP 8.5"]  # the RCP starts from year 2006
Fig_lab = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]# Plot all together
#Fig_lab = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
Magn = ["Low emissions ", "Moderate emissions ", "High emissions "]
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]
   
      
MinYear = 1955
Xlabel = arange(MinYear, 2101., 1.)
            
fsize = 23
fsize2 = 18
Locations_title = ["Great Barrier Reef", "South East Asia", "Caribbean"]

# Parameter resp. for the different region GBR, SEA , CAR
T0_list = array([26.78, 28.13, 27.10]) 
skew_list = array([0.0002, 3.82, 1.06])
rho_list = array([1.0, 0.81, 0.89])
AddTime = 2000*12
count = 1
Col = [(0.827, 0.553, 0.686), (0.416, 0.071, 0.239), (0.69, 0.345, 0.514)]


fsize = 12 #18 #22
fsize2 = 12
MinYear = 1955
count = 1
TempList = linspace(0, 75, 1500)
Gmax = 10./12

fig = plt.figure()
for p in xrange(len(RCP)):
    for z in xrange(len(Locations)):
        sub1 = plt.subplot(3, 3, count)  
        LatC = LatCoordinates[z]
        LonC = LonCoordinates[z]       
        file1 = open("Monthly-SST-scenarios/Months-"+Location[z]+"-"+RCP[p]+"-MPI"+".dat", "r")
        MonthsFinal = load(file1,allow_pickle = True)
        file1.close()
    
        file2 = open("Monthly-SST-scenarios/SST-"+Location[z]+"-"+RCP[p]+"-MPI"+".dat", "r")
        SSTFinal = load(file2,allow_pickle = True)
        file2.close()                
        plt.plot(MonthsFinal, SSTFinal, linewidth = 1, alpha= 0.4, color=(0.7, 0.7, 0.7), label="Used for simulations") 
        T0 = T0_list[z]
        rho = rho_list[z]
        skew = skew_list[z]
        TempCORListCenter = (TempList - T0)/rho
        NormCor = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)
        Dist = Gmax*norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)/max(NormCor)
        maxD = max(Dist)
        T_opt = TempList[int(list(Dist).index(maxD))]
        plt.plot([0, 2150],T_opt*ones(2), "--", linewidth = 1, color = Col[z])                                                                           
        plt.plot(2100+Dist*45, TempList, linewidth = 3,  color=Col[z]) 
        plt.plot(2100*ones(len(TempList)), TempList, linewidth = 0.8,  color="black")                                                 
                                                
        if count not in (7, 8, 9):
            plt.xticks(list(Xlabel[::25])+[2100], [" " for year in list(Xlabel[::25])+[2100]])#, fontsize = fsize2)
        else: 
            plt.xticks(list(Xlabel[::25])+[2100], [int(year) for year in list(Xlabel[::25])+[2100]], rotation=45)#, fontsize = fsize2+2)
            plt.xlabel(r'Year', fontsize = fsize, labelpad=2)
        if count in (1, 4, 7):
            loca = MonthsFinal[0] 
            if RCP_title[p] == "RCP 4.5":
                plt.text(1910, 27, RCP_title[p], rotation="vertical", fontproperties = font, fontsize = fsize)
            else:
                plt.text(1910, 27, RCP_title[p], rotation="vertical", fontproperties = font, fontsize = fsize)
            plt.text(loca-10, 32, Fig_lab[count-1], fontproperties=font, fontsize = fsize2)
            hostTicks = array([0., 0.2, 0.4, 0.6, 0.8, 1.0])
            #plt.yticks(arange(14, 35, 4), fontsize = fsize2)
            plt.ylabel(u"SST (\N{DEGREE SIGN}C)", fontsize = fsize, labelpad = 2)
        elif count in (2, 5, 8):
            plt.text(loca-10, 32, Fig_lab[count-1], fontproperties=font, fontsize = fsize2)
            #plt.yticks(arange(14, 35, 4), [" "]+list([" " for k in arange(14, 35, 4)]), fontsize = fsize2)
        elif count in (3, 6, 9):
            plt.text(loca-10, 31, Fig_lab[count-1], fontproperties=font, fontsize = fsize2)
        if count in (1, 2, 3):
            plt.title(Locations_title[z], fontproperties=font, fontsize = fsize)
        plt.xlim((MinYear, 2100+(45*Gmax)))  
        typ = ":"
        if z == 0:
            plt.ylim((20, 32))
            plt.plot([0, 2150],(T_opt-1.5)*ones(2), typ, linewidth = 1, color = Col[z])                                                                          
            plt.plot([0, 2150],(T_opt+1.5)*ones(2), typ, linewidth = 1, color = Col[z])
            plt.yticks(arange(20, 30+2, 2), ["%d"%sst for sst in arange(20, 30+2, 2)])
            sub1.tick_params(axis = "y", direction = "out", pad = 2)
        elif z == 1:
            plt.ylim((26, 32))
            plt.plot([0, 2150],(T_opt-0.5)*ones(2), typ, linewidth = 1, color = Col[z])                                                                          
            plt.plot([0, 2150],(T_opt+1.5)*ones(2), typ, linewidth = 1, color = Col[z])
            plt.yticks(arange(26, 32, 1), ["%d"%sst for sst in arange(26, 32, 1)])
            sub1.tick_params(axis = "y", direction = "out", pad = 2)

        elif z == 2:
            plt.ylim((25, 31))
            plt.plot([0, 2150],(T_opt-1.5)*ones(2), typ, linewidth = 1, color = Col[z])                                                                          
            plt.plot([0, 2150],(T_opt+1.5)*ones(2), typ, linewidth = 1, color = Col[z])
            plt.yticks(arange(25, 31, 1), ["%d"%sst for sst in arange(25, 31, 1)])
            sub1.tick_params(axis = "y", direction = "out", pad = 2)


        #plt.ylim((18, 32)) 
        # Print warming estimates
        
        compare0 = (MonthsFinal<=2005+11/12)*(MonthsFinal>=1986+11/12)
        compare1 = (MonthsFinal<=2101+11/12)*(MonthsFinal>=2081+11/12)     
        Past = sum(SSTFinal[compare0])/sum(compare0)
        Future = sum(SSTFinal[compare1])/sum(compare1)
        print Locations[z], RCP[p], "Increase in Temperature =", (Future-Past)
          
        count +=1
# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
    
plt.subplots_adjust(bottom = 0.10, right = 0.90, left = 0.20, top = 0.95, wspace = 0.15, hspace = 0.17)      
plt.show()
 
