# -*- coding: utf-8 -*-
from __future__ import division
import pandas as pd
from scipy import arange, mean, array, ones
from matplotlib import pyplot as plt
import numpy as np

# Cover data from Bruno 2007

"""
This code can be used to visualize temporal coral cover data for specific reef location or the regional average coral cover
The Global_Cover_1_5_10.xls data was provided to us by one of the authors (Bruno & Selig 2007) via e-mail.
"""
# Open xls file
import pdb
data0 =   pd.read_excel('Global_Cover_1_5_10.xls', engine='xlrd')
data = data0.copy() 

# Regional location WOD13 temperature data.

GBRlon0, GBRlon1 = (145, 165)
GBRlat0, GBRlat1 = (-28, -10)

SEAlon0, SEAlon1 = (100, 137)
SEAlat0, SEAlat1 = (-10, 13)

CARlon0, CARlon1 = (-80, -65)
CARlat0, CARlat1 = (10, 20)
 

nans = (np.isnan(data.LONGITUDE) & np.isnan(data.LATITUDE))

All_Longs = data.LONGITUDE[~nans]
All_Lats = data.LATITUDE[~nans]
All_Years = data.YEAR[~nans]
All_CoralCover = data.HARD_COR_P[~nans]
All_REEFname = data.REEF_NAME[~nans]


minYear = np.min(All_Years)
maxYear = np.max(All_Years)

GBR_Years = []
SEA_Years = []
CAR_Years = []

GBR_Cover = []
SEA_Cover = []
CAR_Cover = []

GBR_Reefs = []
SEA_Reefs = []
CAR_Reefs = []

GBR_Cover_Sep = []
SEA_Cover_Sep = []
CAR_Cover_Sep = []

GBR_Years_Sep = []
SEA_Years_Sep = []
CAR_Years_Sep = []


Year_list = arange(minYear, maxYear + 1, 1., dtype = "int64") # int64 is the type of the numbers in data.YEAR

for i in range(len(Year_list)):
    IndexGBR_i = ((All_Longs >= GBRlon0)&(All_Longs < GBRlon1)&(All_Lats >= GBRlat0)&(All_Lats < GBRlat1)&(All_Years == Year_list[i]))
    IndexSEA_i = ((All_Longs >= SEAlon0)&(All_Longs < SEAlon1)&(All_Lats >= SEAlat0)&(All_Lats < SEAlat1)&(All_Years == Year_list[i]))
    IndexCAR_i = ((All_Longs >= CARlon0)&(All_Longs < CARlon1)&(All_Lats >= CARlat0)&(All_Lats < CARlat1)&(All_Years == Year_list[i]))
    
    if sum(IndexGBR_i) !=0: 
        GBR_Years.append(Year_list[i])
        All_Cover = All_CoralCover[IndexGBR_i]
        # Include unique measurement per reef name, assuming the % cover reflect the whole reef area and measurements could be made by different people on different months of the same year 
        Names = All_REEFname[IndexGBR_i] 
        Names = np.array([n.encode('utf-8') for n in Names])
        UniqueReefs = set(list(Names)) 
        UniqueMes = []
        for reef in UniqueReefs:
            allMes = All_Cover[Names == reef] # take all measurement at a particular reef
            UniqueMes.append(mean(allMes))    
        GBR_Cover.append(mean(UniqueMes)) 
        
        GBR_Reefs = GBR_Reefs + list(UniqueReefs)
        GBR_Cover_Sep = GBR_Cover_Sep + UniqueMes
        GBR_Years_Sep = GBR_Years_Sep +   list(Year_list[i]*ones(len(UniqueMes)))
        #plt.plot(Year_list[i]*ones(len(UniqueMes)), array(UniqueMes), "o")
        # Ignore multiple measurements per reef name, assuming the % cover does not reflect the whole reef area but only some subset of it, there is not reference to surface area of the reef in the dataset
        #GBR_Cover.append(mean(All_CoralCover[IndexGBR_i[0]]))

    if sum(IndexSEA_i) !=0:
        SEA_Years.append(Year_list[i])
        All_Cover = All_CoralCover[IndexSEA_i]
        Names = All_REEFname[IndexSEA_i]
        Names = np.array([n.encode('utf-8') for n in Names])
        UniqueReefs = set(list(Names))
        UniqueMes = []
        for reef in UniqueReefs:
            try:
                allMes = All_Cover[Names == reef]
            except:
                pdb.set_trace()
            UniqueMes.append(mean(allMes))
        SEA_Cover.append(mean(UniqueMes))
        SEA_Reefs = SEA_Reefs + list(UniqueReefs)
        SEA_Cover_Sep = SEA_Cover_Sep + UniqueMes
        SEA_Years_Sep = SEA_Years_Sep + list(Year_list[i]*ones(len(UniqueMes)))
        #SEA_Cover.append(mean(All_CoralCover[IndexSEA_i[0]]))
    if sum(IndexCAR_i) !=0:
        CAR_Years.append(Year_list[i])
        All_Cover = All_CoralCover[IndexCAR_i]
        Names = All_REEFname[IndexCAR_i]
        Names = np.array([n.encode('utf-8') for n in Names])
        UniqueReefs = set(list(Names))
        UniqueMes = []
        for reef in UniqueReefs:
            allMes = All_Cover[Names == reef]
            UniqueMes.append(mean(allMes))
        CAR_Cover.append(mean(UniqueMes))
        CAR_Reefs = CAR_Reefs + list(UniqueReefs)
        CAR_Cover_Sep = CAR_Cover_Sep + UniqueMes
        CAR_Years_Sep = CAR_Years_Sep +   list(Year_list[i]*ones(len(UniqueMes)))
        #CAR_Cover.append(mean(All_CoralCover[IndexCAR_i[0]]))


# Plot yearly dataset for specific reefs 
Chosen_set_GBR = ["Lizard Island", "Davies", "Low Isles", "Chinaman", "Turner Cay", "Linnet", "Hayman Island", "Heron Island"]
Chosen_set_SEA = ["South Bang Tao", "North Bang Tao", "North Patong", "Bantayan Beach", "Campuyo Reef"]
Chosen_set_CAR = ["West Flower Garden Banks", "Chengue Bay", "west forereef", "Yawzi Point", "Newfound Reef"]
REG_Chosen_set = [Chosen_set_GBR, Chosen_set_SEA, Chosen_set_CAR]

REG = ["GBR", "SEA", "CAR"]

REG_Cover = [GBR_Cover, SEA_Cover, CAR_Cover]
REG_Reefs = [GBR_Reefs, SEA_Reefs, CAR_Reefs]
REG_Years = [GBR_Years, SEA_Years, CAR_Years]
REG_Cover_Sep = [GBR_Cover_Sep, SEA_Cover_Sep, CAR_Cover_Sep]
REG_Years_Sep = [GBR_Years_Sep, SEA_Years_Sep, CAR_Years_Sep]

MinData = [10, 6, 6]
s = 1
plt.figure(figsize=(90, 60))

for k in range(3):
    Unique_Reg_Reefs = set(list(REG_Reefs[k]))   
    
    Reg_Reefs = array(REG_Reefs[k])
    Reg_Cover_Sep = array(REG_Cover_Sep[k])
    Reg_Years_Sep = array(REG_Years_Sep[k])  
    
    Chosen_set = REG_Chosen_set[k]
    Reg_Years = REG_Years[k]
    
    Target_reefs = []
    plt.subplot(3,2, s)
    num = 0
    lim = 0
    for reef in Unique_Reg_Reefs:
        indexReef_i = (Reg_Reefs == reef).nonzero()
        if len(indexReef_i[0])!=0: # the first element is what we need
            if len(indexReef_i[0])>=MinData[k] :#len(GBR_Years)/2: # choose chosen_set from reefs with several years measurments
                #print reef  # in GBR 277 reefs measured less than 5 times,46 reefs measured more than or =10 years, 42 measured more than or = 12 years
                Target_reefs.append(reef)
                num +=1   
            ReefCov = Reg_Cover_Sep[indexReef_i[0]]
            YearsReefCov = Reg_Years_Sep[indexReef_i[0]]
            if reef in Chosen_set and len(indexReef_i[0])>=MinData[k] and lim < 4:
                lim += 1 # limit the number of reef highlighted
                plt.plot(YearsReefCov, ReefCov, "o", markersize = 5, alpha = 0.5, color = "black")
                plt.plot(YearsReefCov, ReefCov, linewidth = 2, label = reef)#, color = (0.063, 0.302, 0.))
            elif len(indexReef_i[0])>=MinData[k]:
                plt.plot(YearsReefCov, ReefCov, "o", markersize = 5, alpha = 0.1, color = "grey")
                plt.plot(YearsReefCov, ReefCov, linewidth = 0.1, color = "black")               
    print("Num"+REG[k], num)
    print(Target_reefs)
    plt.title("Target reef in the "+REG[k]+" (measured $\geq$ %d times, %d reefs out of %d)"%(MinData[k], num, len(Unique_Reg_Reefs)), fontsize = 10)
    plt.xlim((1970, 2007))
    if s == 5:
        plt.xticks(arange(1970, 2007, 5))
    else:
        plt.xticks(arange(1970, 2007, 5), ["" for a in arange(1970, 2007, 5)])
    plt.ylim(0, 100)      

    plt.legend(loc="upper left", fontsize = 7)
    
    
    # Take the mean coral cover accross fixed target reefs 
    Reg_Cover_target = []
    Year_target = [] # year with all target reefs measured
    for i in range(len(Reg_Years)):
        Cov_mes = []
        for reef in Target_reefs:
            indexReef_i = (Reg_Reefs == reef).nonzero()
            ReefCov = Reg_Cover_Sep[indexReef_i[0]]
            YearsReefCov = Reg_Years_Sep[indexReef_i[0]]
            if Reg_Years[i] in YearsReefCov: # take the cover measurment for that year
                Cov_mes.append(ReefCov[list(YearsReefCov).index(Reg_Years[i])]) 
            #else: # assume that the cover is the mean of all existing measurement for that reef
            #Cov_mes.append(mean(ReefCov))
        if len(Cov_mes) == len(Target_reefs):
            Reg_Cover_target.append(mean(Cov_mes))
            Year_target.append(Reg_Years[i])
    
    plt.ylabel("Coral cover ($\%$)", fontsize = 10) 
    plt.xlabel("Year", fontsize = 10)       
    plt.subplot(3, 2, s+1)  
    plt.title("Mean coral cover: all target reef (red), all reefs (green)", fontsize = 10)        
    plt.plot(Year_target, Reg_Cover_target, "o", markersize = 5, alpha = 0.5, color = "red", label="All target reefs measured")
    plt.xlim((1970, 2007))
    if s == 5:
        plt.xticks(arange(1970, 2007, 5))
    else:
        plt.xticks(arange(1970, 2007, 5), ["" for a in arange(1970, 2007, 5)])     
    
    Reg_Cover = REG_Cover[k]
    plt.plot(Reg_Years, Reg_Cover, "o", markersize = 5, alpha = 0.5, color = (0.063, 0.302, 0.), label="Mean of all reef meas.")
    plt.xticks(arange(1970, 2007, 5))
    plt.xlabel("Year", fontsize = 15) 
    plt.ylim(0, 100)      
    
    s +=2
 
             
                           
# Output for Regional Average Coral Cover
"""
#GBR
GBR_Years = [1971, 1972, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 
1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007]


# Include unique measurement per reef name
GBR_Cover_1 = [4.2999999999999998, 34.975000000000001, 49.375, 24.906481481481478, 24.441163003663007, 
24.126298701298698, 16.893939393939394, 22.58467261904762,17.311111111111114, 21.028985714285714, 
23.527606393606391, 25.455069209922868, 28.285919266492023, 31.795480366220563, 32.272448584582911,
32.982510303536621, 33.638900123685843, 29.743189587198515, 30.616005555555553, 27.787401217829789,
31.235181694495974, 30.107761904761908, 26.650676470588238, 30.476910112359551,14.855399999999999]

# Ignore multiple measurements per reef name
GBR_Cover_2 = [4.2999999999999998, 34.975000000000001, 49.375, 25.204545454545453, 23.355555555555554,
25.05263157894737, 14.851851851851851, 22.627118644067796, 15.0, 20.949863636363638, 23.404128205128206, 
25.099387499999999, 26.204821428571432, 28.828654929577468, 29.184980519480519, 29.911180327868852, 
30.113549668874178, 27.198119402985078, 27.90510344827586, 26.422128834355828, 28.465675862068966, 
29.903699029126212, 26.18496, 30.869656250000002, 14.855399999999999]


plt.subplot(1, 3, 1)
plt.title("Great Barrier Reef")
plt.plot(GBR_Years, GBR_Cover_1, "o", markersize = 5, alpha = 0.5, color = (0.063, 0.302, 0.), label = "unique mes per reef name")
plt.plot(GBR_Years, GBR_Cover_2, "o", markersize = 5, alpha = 0.5, color = "purple", label = "multiple mes per reef name")

plt.plot(GBR_Years, GBR_Cover_1, linewidth = 1, markersize = 5, alpha = 1, color = (0.063, 0.302, 0.))
plt.plot(GBR_Years, GBR_Cover_2, linewidth = 1, markersize = 5, alpha = 1, color = "purple")

plt.ylim((0, 60))


#SEA
SEA_Years = [1971, 1974, 1977, 1978, 1979, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 
1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006]

# Include unique measurement per reef name
SEA_Cover_1 = [65.0, 39.683333333333337, 15.176470588235293, 10.75, 40.615229885057467, 38.426666666666669, 47.682051282051283,
29.899999999999999, 34.222222222222221, 27.421666666666663, 38.978494623655912, 42.787939221272559, 35.265217391304347,
30.222222222222221, 49.857142857142854, 38.984946236559139, 30.648148148148149, 35.902311594202892, 39.000386904761903, 
38.747898397435897, 27.335411368015414, 34.826259328358212, 27.322793103448276, 29.019946808510639, 30.414325842696631,
33.062789351851855, 31.652665770609318, 34.652255639097746, 31.402439024390244, 27.543560606060606, 33.609910337552741]

# Ignore multiple measurement per reef name
SEA_Cover_2 = [65.0, 39.6875, 15.176470588235293, 10.75, 39.769148936170218, 37.195, 43.788888888888891, 30.3125, 
34.222222222222221, 25.473529411764705, 38.348837209302324, 36.396551724137929, 37.450000000000003, 33.272727272727273, 
51.291666666666664, 38.555555555555557, 29.657894736842106, 36.025822916666669, 37.633065217391298, 41.574519148936162, 
30.490393023255812, 35.980434782608697, 29.448816720257234, 28.265822784810126, 29.483267716535433, 33.163580246913583, 
30.846153846153847, 33.326086956521742, 31.410845588235293, 26.388888888888889, 33.63180434782609]

plt.subplot(1, 3, 2)
plt.title("South East Asia")
plt.plot(SEA_Years, SEA_Cover_1, "o", markersize = 5, alpha = 0.5, color = (0.063, 0.302, 0.), label = "unique mes per reef name")
plt.plot(SEA_Years, SEA_Cover_2, "o", markersize = 5, alpha = 0.5, color = "purple", label = "multiple mes per reef name")

plt.plot(SEA_Years, SEA_Cover_1, linewidth = 1, alpha = 1, color = (0.063, 0.302, 0.))
plt.plot(SEA_Years, SEA_Cover_2, linewidth = 1, alpha = 1, color = "purple")
plt.ylim((0, 60))

#CAR
CAR_Years = [1976, 1977, 1978, 1980, 1981, 1982, 1983, 1984, 1986, 1987, 1988, 1989,1990, 1992,
1993, 1994, 1995,1996,1997,1998,1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006]

# Include unique measurement per reef name
CAR_Cover_1 = [25.17, 47.1, 46.7735, 28.0875, 16.1, 22.4429375,  17.2,  12.26666667,
         6.78571429, 10.296, 22.4, 2.55, 27.85,   6.3375,  22.32375, 24.10105556,
        26.1696, 17.36739394, 25.76126312, 32.40814444, 28.30322, 15.62805952, 20.24744444, 35.37140909,
        26.1427, 24.74643878,21.81088542, 14.9610625 ]

# Ignore multiple measurement per reef name
CAR_Cover_2 = [25.170000000000002, 47.100000000000001, 39.437600000000003, 34.020000000000003, 16.100000000000001, 
19.676333333333332,17.199999999999999, 12.266666666666666, 6.7857142857142856, 10.295999999999999,22.399999999999999,
2.5499999999999998, 27.850000000000001, 6.3374999999999995, 22.323749999999997, 19.444875, 23.513714285714286,
11.740631578947369,18.934854545454549, 32.323315789473689, 28.57958620689655,17.889960784313725,19.355499999999999,
29.202666666666669, 25.908627118644063, 24.661862068965515, 20.324777777777779, 15.035613636363637]

plt.subplot(1, 3, 3)
plt.title("Caribbean")
plt.plot(CAR_Years, CAR_Cover_1, "o", markersize = 5, alpha = 0.5, color = (0.063, 0.302, 0.), label = "unique mes per reef name")
plt.plot(CAR_Years, CAR_Cover_2, "o", markersize = 5, alpha = 0.5, color = "purple", label = "multiple mes per reef name")

plt.plot(CAR_Years, CAR_Cover_1, linewidth = 1, alpha = 1, color = (0.063, 0.302, 0.))
plt.plot(CAR_Years, CAR_Cover_2, linewidth = 1, alpha = 1, color = "purple")
plt.ylim((0, 60))

plt.legend(loc = "upper right", ncol = 1)
"""
# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
    
    