# -*- coding: utf-8 -*-
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from scipy import linspace
import matplotlib as mpl

import pandas as pd
from scipy import arange, mean, array, ones, zeros, concatenate, std
from matplotlib import pyplot as plt
from scipy.stats import norm
import pdb

"""
This code was used to generate Coral Cover Data in manuscript: Appendix Figure 2
The Global_Cover_1_5_10.xls data was provided to us by Bruno 2017 via e-mail.
"""
# Open xls file
data0 =   pd.read_excel('Global_Cover_1_5_10.xls', engine='xlrd')
data = data0.copy() 


from matplotlib.font_manager import FontProperties
font = {'family' : 'normal',
        'weight' : 'bold'}

fsize = 16#20
fsize2 = 13#13
fsize3 = 10
plt.rcParams["font.family"] = "arial"   

font0 = FontProperties()
font = font0.copy()
weight = "bold"
font.set_weight(weight)    
                                                                                                                                                                                                                                                                                                                                                     
fig = plt.figure()

Col = [(0.827, 0.553, 0.686), (0.69, 0.345, 0.514), (0.416, 0.071, 0.239)]
#landColor = (0.286, 0.686, 0.408)
landColor = (0.7, 0.7, 0.7) #(0.28, 0.68, 0.402)
#oceanColor = (0.55, 0.897, 0.975)#(0.62, 0.745, 0.839)
oceanColor = (1, 1, 1)#(0.62, 0.745, 0.839)

# GBR (Gread Barrier Reef)
sub1 = fig.add_subplot(4, 3, 1)  
plt.title("Great Barrier Reef", fontsize = fsize, fontproperties = font)
                                                                                                                                                                              
my_map = Basemap(projection='cyl')  # longitudes are always between -pi, pi, lat = -pi/2, pi/2
#my_map.drawcoastlines()
#my_map.drawmapboundary(color = (0.427, 0.573, 0.627)) # deprecated
sub1.set_facecolor(oceanColor)
my_map.fillcontinents(color = landColor)

GBRx1,GBRy1 = my_map(145, -28) 
GBRx2,GBRy2 = my_map(145, -10)
GBRx3,GBRy3 = my_map(165, -10)
GBRx4,GBRy4 = my_map(165, -28)

# Just GBR close to coastline
plt.xticks(arange(145, 165+5, 5), [u"%d\N{DEGREE SIGN}E"%ar for ar in arange(145, 165+5, 5)], fontsize = fsize3)
plt.yticks(arange(-28, -10+2, 2), [u"%d\N{DEGREE SIGN}S"%abs(ar) for ar in arange(-28, -10+2, 2)], fontsize = fsize3)

plt.xlim((142, 167))
plt.ylim((-28, -12))

# Just GBR zoom out
#plt.xlim(120, 170)
#plt.ylim(-35, -5)

#SEA‚ Southeast Asia
sub2 = fig.add_subplot(4, 3, 2)
plt.title("South East Asia", fontsize = fsize, fontproperties = font)
                                                                                                                                                                              
my_map = Basemap(projection='cyl')  # longitudes are always between -pi, pi, lat = -pi/2, pi/2
#my_map.drawcoastlines()
sub2.set_facecolor(oceanColor)
my_map.fillcontinents(color = landColor)

SEAx1,SEAy1 = my_map(100, -10)
SEAx2,SEAy2 = my_map(100, 13)
SEAx3,SEAy3 = my_map(137, 13)
SEAx4,SEAy4 = my_map(137, -10)

# Just SEA
plt.xticks(arange(100, 140+10, 10), [u"%d\N{DEGREE SIGN}E"%ar for ar in arange(100, 140+10, 10)], fontsize = fsize3)
plt.yticks(list(arange(-12, 0, 4))+list(arange(0, 16+4, 4)), [u"%d\N{DEGREE SIGN}S"%abs(ar) for ar in arange(-12, 0, 4)]+[u"%d\N{DEGREE SIGN}N"%abs(ar) for ar in arange(0, 16+4, 4)], fontsize = fsize3)
plt.xlim((95, 142))
plt.ylim((-12, 16))

# CAR (Caribbean) 
sub3 = fig.add_subplot(4, 3, 3) 
plt.title("Caribbean", fontsize = fsize, fontproperties = font)                                                                                                                                                                            
my_map = Basemap(projection='cyl')  # longitudes are always between -pi, pi, lat = -pi/2, pi/2
#my_map.drawcoastlines()
sub3.set_facecolor(oceanColor)
my_map.fillcontinents(color = landColor)
CARx1,CARy1 = my_map(-80, 10)    
CARx2,CARy2 = my_map(-80, 20)
CARx3,CARy3 = my_map(-65, 20)
CARx4,CARy4 = my_map(-65, 10)

# Just CAR
plt.xticks(arange(-80, -65+5, 5), [u"%d\N{DEGREE SIGN}W"%abs(ar) for ar in arange(-80,65+5, 5)], fontsize = fsize3)
plt.yticks(arange(8, 22+2, 2), [u"%d\N{DEGREE SIGN}N"%abs(ar) for ar in arange(8, 22+2, 2)], fontsize = fsize3)

plt.xlim((-82, -62))
plt.ylim((8, 22))


# get the relevant data
nans = (np.isnan(data.LONGITUDE) & np.isnan(data.LATITUDE))
All_Longs = data.LONGITUDE[~nans]
All_Lats = data.LATITUDE[~nans]
All_Years = data.YEAR[~nans]
All_CoralCover = data.HARD_COR_P[~nans]
All_REEFname = data.REEF_NAME[~nans]


GBRlon0, GBRlon1 = (145, 165)
GBRlat0, GBRlat1 = (-28, -10)

SEAlon0, SEAlon1 = (100, 137)
SEAlat0, SEAlat1 = (-10, 13)

CARlon0, CARlon1 = (-80, -65)
CARlat0, CARlat1 = (10, 20)

minYear = min(All_Years)
maxYear = max(All_Years)

GBR_Years = []
SEA_Years = []
CAR_Years = []

GBR_Infos = []
SEA_Infos = []
CAR_Infos = []


Year_list = arange(minYear, maxYear + 1, 1., dtype = "int64") # int64 is the type of the numbers in data.YEAR

source2Col = (0.353, 0.675, 0.337)
for i in xrange(len(Year_list)):
    IndexGBR_i = ((All_Longs >= GBRlon0)&(All_Longs < GBRlon1)&(All_Lats >= GBRlat0)&(All_Lats < GBRlat1)&(All_Years == Year_list[i]))
    IndexSEA_i = ((All_Longs >= SEAlon0)&(All_Longs < SEAlon1)&(All_Lats >= SEAlat0)&(All_Lats < SEAlat1)&(All_Years == Year_list[i]))
    IndexCAR_i = ((All_Longs >= CARlon0)&(All_Longs < CARlon1)&(All_Lats >= CARlat0)&(All_Lats < CARlat1)&(All_Years == Year_list[i]))
    if sum((IndexGBR_i)) != 0: 
        GBR_Years.append(Year_list[i])
        GBR_Infos.append([All_Longs[IndexGBR_i], All_Lats[IndexGBR_i], All_CoralCover[IndexGBR_i], All_REEFname[IndexGBR_i]])
    if sum((IndexSEA_i)) !=0:
        SEA_Years.append(Year_list[i])
        SEA_Infos.append([All_Longs[IndexSEA_i], All_Lats[IndexSEA_i], All_CoralCover[IndexSEA_i], All_REEFname[IndexSEA_i]])
    if sum((IndexCAR_i)) !=0:
        CAR_Years.append(Year_list[i])
        CAR_Infos.append([All_Longs[IndexCAR_i], All_Lats[IndexCAR_i], All_CoralCover[IndexCAR_i], All_REEFname[IndexCAR_i]])

# Coral cover Vs Year 
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# GBR
sub4 = fig.add_subplot(4, 3, 4)
datm = 3 # data markersize

NumbReefGBR = zeros(len(GBR_Years))
FarReef = []
for m in xrange(len(GBR_Years)):
    partGBR = GBR_Infos[m]
    Longs = array(partGBR[0])
    Lats = array(partGBR[1])
    Cover = array(partGBR[2])
    Names = array(partGBR[3])
    Names = np.array([n.encode('utf-8') for n in Names])
    UniqueReef = set(list(Names)) # deals with repeated measurments of different samples on a reef
    NumbReefGBR[m] = len(UniqueReef)
    UniqueMes = []
    UniqueLong = []
    UniqueLat = []
    for reef in UniqueReef:
        allMes = Cover[Names == reef] # take all measurment at a particular reef
        UniqueMes.append(mean(allMes))
        UniqueLong.append(list(Longs[Names == reef])[0]) # just take the first 'cause they are all the same for each unique reef, I already checked
        UniqueLat.append(list(Lats[Names == reef])[0])
        if list(Longs[Names == reef])[0]>153 and reef not in FarReef: # limits in Death et al. 2012
            FarReef.append(reef)
            print reef, list(Longs[Names == reef])[0], list(Lats[Names == reef])[0]
            
    x_coord = ones(len(UniqueReef))*GBR_Years[m]
    plt.plot(x_coord, UniqueMes, "o", markersize = datm, color = source2Col, alpha = 0.3)
    plt.plot(x_coord, UniqueMes, "-", markersize = datm, linewidth = 0.25, color = source2Col, alpha = 0.3)
    # Measurement location for all years
    
    pLongs = list(UniqueLong)
    pLats = list(UniqueLat)
    pCover = list(UniqueMes)
    for k in xrange(len(pLongs)):
        sub1.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    
    """
    # Measurement location in a particular year
    if GBR_Years[m] ==1987:
        pLongs = list(UniqueLong)
        pLats = list(UniqueLat)
        pCover = list(UniqueMes)
        for k in xrange(len(pLongs)):
            sub1.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
        sub4.plot(x_coord, UniqueMes, "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    """    
sub4.set_ylabel("Coral Cover", color = source2Col, fontsize = fsize)
sub4.set_yticks(((0, 20, 40, 60, 80, 100)))
sub4.set_yticklabels((0, 20, 40, 60, 80, 100), color = source2Col, fontsize = fsize2)
sub4.set_ylim((0, 100))
sub4.set_xticks(arange(1970, 2020, 10))
sub4.set_xticklabels([" "%k for k in arange(1970, 2020, 10)], fontsize = fsize2)
sub4.set_xlim ((1970, 2010))
plt.xlabel("Year", fontsize = fsize)

plt1 = sub4.twinx()
plt1.plot(GBR_Years, NumbReefGBR, linewidth = 2, color = (0.502, 0.082, 0.082), alpha = 0.75)
plt1.set_yticks(arange(0, 120, 20))
plt1.set_yticklabels(arange(0, 120, 20), color = (0.502, 0.082, 0.082), fontsize = fsize2)
plt1.set_ylim((0, 100))
plt1.set_xticks(arange(1970, 2020, 10))
plt1.set_xticklabels(["%d"%k for k in arange(1970, 2020, 10)], fontsize = 9)
#plt1.set_xticklabels([" "%k for k in arange(1970, 2020, 5)], fontsize = fsize2)
plt1.set_xlim ((1970, 2010))

# SEA                       
sub5 = fig.add_subplot(4, 3, 5)
plt.yticks((0, 20, 40, 60, 80, 100), ("", "", "", "", "", "", ""), fontsize = fsize)
plt.xticks(arange(1970, 2015, 5), [" "%k for k in arange(1970, 2015,5)], fontsize = fsize2)
plt.xlim ((1970, 2010))
plt.ylim((0, 100))
plt.xlabel("Year", fontsize = fsize)
NumbReefSEA = zeros(len(SEA_Years))
for m in xrange(len(SEA_Years)):
    partSEA = SEA_Infos[m]
    Longs = array(partSEA[0])
    Lats = array(partSEA[1])
    Cover = array(partSEA[2])
    Names = array(partSEA[3])
    Names = np.array([n.encode('utf-8') for n in Names])
    UniqueReef = set(list(Names)) # deals with repeated measurments
    NumbReefSEA[m] = len(UniqueReef)
    UniqueMes = []
    UniqueLong = []
    UniqueLat = []
    for reef in UniqueReef:
        allMes = Cover[Names == reef] # take all measurment at a particular reef
        UniqueMes.append(mean(allMes))
        UniqueLong.append(list(Longs[Names == reef])[0]) # just take the first 'cause they are all the same for each unique reef
        UniqueLat.append(list(Lats[Names == reef])[0])
    x_coord = ones(len(UniqueReef))*SEA_Years[m]
    plt.plot(x_coord, UniqueMes, "o", markersize = datm, color = source2Col, alpha = 0.3)
    plt.plot(x_coord, UniqueMes, "-", markersize = datm, linewidth = 0.25, color = source2Col, alpha = 0.3)
    # Measurement location for all years
    pLongs = list(UniqueLong)
    pLats = list(UniqueLat)
    pCover = list(UniqueMes)
    for k in xrange(len(pLongs)):
        sub2.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    
    """
    # Measurement locations in a particular year 
    if SEA_Years[m] == 2000:
        pLongs = list(UniqueLong)
        pLats = list(UniqueLat)
        pCover = list(UniqueMes)
        for k in xrange(len(pLongs)):
            sub2.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
        # Associated coral cover data
        sub5.plot(x_coord, UniqueMes, "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    """
plt2 = sub5.twinx()
plt2.plot(SEA_Years, NumbReefSEA, linewidth = 2, color = (0.502, 0.082, 0.082), alpha = 0.75)
plt2.set_yticks(arange(0, 400, 50))
plt2.set_yticklabels(arange(0, 400, 50), color = (0.502, 0.082, 0.082), fontsize = fsize2)
plt2.set_ylim((0, 350))
plt2.set_xticks(arange(1970, 2020, 10))
plt2.set_xticklabels(["%d"%k for k in arange(1970, 2020, 10)], fontsize = 9)
#plt2.set_xticklabels([" "%k for k in arange(1970, 2020, 5)], fontsize = fsize2)
plt2.set_xlim ((1970, 2010))

# CAR
sub6 = fig.add_subplot(4, 3, 6)
plt.yticks((0, 20, 40, 60, 80, 100), ("", "", "", "", "", "", ""), fontsize = fsize)
plt.xticks(arange(1970, 2015, 10), [" "%k for k in arange(1970, 2015, 10)], fontsize = fsize2)
plt.xlim ((1970, 2010))
plt.ylim((0, 100))
plt.xlabel("Year", fontsize = fsize)

NumbReefCAR = zeros(len(CAR_Years))
for m in xrange(len(CAR_Years)):
    partCAR = CAR_Infos[m]
    Longs = array(partCAR[0])
    Lats = array(partCAR[1])
    Cover = array(partCAR[2])
    Names = array(partCAR[3])
    Names = np.array([n.encode('utf-8') for n in Names])
    UniqueReef = set(list(Names)) # deals with repeated measurments
    NumbReefCAR[m] = len(UniqueReef)
    UniqueMes = []
    UniqueLong = []
    UniqueLat = []
    for reef in UniqueReef:
        allMes = Cover[Names == reef] # take all measurment at a particular reef
        UniqueMes.append(mean(allMes))
        UniqueLong.append(list(Longs[Names == reef])[0]) # just take the first 'cause they are all the same for each unique reef
        UniqueLat.append(list(Lats[Names == reef])[0])
    x_coord = ones(len(UniqueReef))*CAR_Years[m]
    plt.plot(x_coord, UniqueMes, "o", markersize = datm, color = source2Col, alpha = 0.3)
    plt.plot(x_coord, UniqueMes, "-", markersize = datm, linewidth = 0.25, color = source2Col, alpha = 0.3)
    # Measurement location for all years
    pLongs = list(UniqueLong)
    pLats = list(UniqueLat)
    pCover = list(UniqueMes)
    for k in xrange(len(pLongs)):
        sub3.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    
    """
    # Measurment location in a particular year 
    if CAR_Years[m] == 2000:
        pLongs = list(UniqueLong)
        pLats = list(UniqueLat)
        pCover = list(UniqueMes)
        for k in xrange(len(pLongs)):
            sub3.plot(pLongs[k], pLats[k], "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
        # Associated coral cover data
        sub6.plot(x_coord, UniqueMes, "o", color = source2Col, markersize = 4, markeredgecolor = "black", markeredgewidth = 0.5)
    """
plt3 = sub6.twinx()
plt3.plot(CAR_Years, NumbReefCAR, linewidth = 2, color = (0.502, 0.082, 0.082), alpha = 0.75)
plt3.set_yticks(arange(0, 120, 20))
plt3.set_yticklabels(arange(0, 120, 20), color = (0.502, 0.082, 0.082), fontsize = fsize2)
plt3.set_ylim((0, 100))
plt3.set_ylabel("Number of \n Reef Surveyed", color = (0.502, 0.082, 0.082), fontsize = fsize)
plt3.set_xticks(arange(1970, 2020, 10))
plt3.set_xticklabels(["%d"%k for k in arange(1970, 2020, 10)], fontsize = 9)

#plt3.set_xticklabels([" "%k for k in arange(1970, 2020, 5)], fontsize = fsize2)
plt3.set_xlim ((1970, 2010))


# Plot with maximal window
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()    
plt.show()