from __future__ import division
from matplotlib import pyplot as plt
from scipy import array, linspace, ones, concatenate, mean, arange
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
plt.switch_backend('agg')

""" 
Code producing Symbiont Biomass reduction by Coral Type and by Location (Appendix Figure S1)

"""

# T for temperature increase and S for associated symbiont densities as % compared from healthy or unbleached corals
# Reference: Table: S-density-bleaching
# Ti_initial could be monthly/yearly average temperature on the site or in control experiment for unbleached corals column 4 of Table: S-density-bleaching
# Ti_MMMax is the value we could asign the the mean monthly maximum temperature used as reference for the conventional DHM or DHW. We chose this value to be 
# still within the range of reported measurements and/or for which there is associated symbiont density (column 4 of Table: S-density-bleaching)

plt.rcParams["font.family"] = "arial"         
                                                                                                                                                                                                                                                                                                                                                     
#plt.subplots_adjust(hspace=0.0001)                                                                                                                                                                                                                                                
#plt.subplot(2, 1, 2)

fsize = 18*2.5#20
fsize2 = 16*2.5#13
plt.rcParams["font.family"] = "arial"   

font0 = FontProperties()
font = font0.copy()
weight = "bold"
font.set_weight(weight) 

fig = plt.figure(figsize = (35, 15))
sub2 = fig.add_subplot(1, 2, 1)
plt.title("By coral types\n", fontsize = fsize, fontproperties = font)
plt.ylabel("Reduction in symbiont biomass \n $\Delta S$ (%)", fontsize = fsize)
plt.xlabel("$\Delta T$"+u" (\N{DEGREE SIGN}C)", fontsize = fsize)
plt.xticks((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), fontsize = fsize2)
plt.text(-1.2, 101, "a", fontproperties = font, fontsize = fsize2)

plt.xlim(0, 9)

sub2.tick_params(axis = "y", direction = "in", labelsize = fsize2)
#sub2.set_yticks((0, 25, 50, 75, 100)) # GBR and SEA
ylabs = arange(0, 110, 10)
sub2.set_yticks(ylabs) # GBR and SEA
sub2.set_yticklabels(list(["%d"%s for s in ylabs]))
sub2.set_ylim((0, 100))
part2 = sub2.twinx()
part2.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
part2.set_yticks(ylabs)
part2.set_yticklabels([" " for d in ylabs])
part2.set_ylim((0, 100))


sub3 = fig.add_subplot(1, 2, 2)
plt.title("By study regions\n", fontsize = fsize, fontproperties = font)
plt.xlabel("$\Delta T$"+u" (\N{DEGREE SIGN}C)", fontsize = fsize)
plt.xticks((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), fontsize = fsize2)
plt.xlim(0, 9)
plt.text(-0.5, 101, "b", fontproperties = font, fontsize = fsize2)

sub3.tick_params(axis = "y", direction = "in", labelsize = fsize2)
sub3.set_yticks(ylabs) # GBR and SEA
sub3.set_yticklabels(list([" "%s for s in ylabs]))
sub3.set_ylim((0, 100))
part3 = sub3.twinx()
part3.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
part3.set_yticks(ylabs)
part3.set_yticklabels([" " for d in ylabs])
part3.set_ylim((0, 100))

a0 = 1
m0 = 20

#color for references
numLoc = 16
cmx = mpl.cm
cmapC = mpl.cm.tab20c#Vega20c
cNorm = mpl.colors.Normalize(vmin=1, vmax=numLoc)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmapC)

colorbranching = "blue"
colormassive = "orange"
colorotherShape = "grey"

Loc_target = ["GBR", "SEA", "CAR"]
#col = [(0.827, 0.553, 0.686), (0.416, 0.071, 0.239), (0.69, 0.345, 0.514), "grey"]
col = ["green",(0.44, 0.07, 0.20), (0.8, 0.4, 0.8), (0.9, 0.9, 0.9)]  
label_unique = [] # for regions

#[1] separate between those with visual signs and those with no visual signs of bleaching
S1 = array([(10+11)/2, (28+16)/2, (52+47)/2, (55+69)/2])
T1 = (30.8)*ones(len(S1))
T1_initial = 28.6  
T1_MMMax = 29.8

S1branching = array([(10+11)/2, (52+47)/2])
T1branching = (30.8)*ones(len(S1branching))

S1massive = array([(28+16)/2, (55+69)/2])
T1massive = (30.8)*ones(len(S1massive))

S1Location = "Other"#"South China Sea"

#sub2.plot(T1branching-T1_initial, S1branching, "o", markersize = m0, alpha = a0, color = colorbranching, label = "Branching")

#sub2.plot(T1massive-T1_initial, S1massive, "o", markersize = m0, alpha = a0, color = colormassive, label = "Massive")

if S1Location not in label_unique and S1Location in Loc_target:
    sub3.plot(T1-T1_initial, 100-S1, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S1Location)], label = S1Location)
    label_unique.append(S1Location)
elif S1Location in Loc_target:
    sub3.plot(T1-T1_initial, 100-S1, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S1Location)])
    
#[2]
S2 = 100*array([1.5/3.5, 0.9/1.5, 2.5/5, 1/3.5, 2/5, 1/3, 1/1.5, 2/2.5, 1.5/2])
T2 = ((30+31+32)/3) * ones(len(S2))
T2_initial = (26 + 27)/2
T2_MMMax = 29

S2branching = 100*array([1.5/3.5, 0.9/1.5, 1/1.5])
T2branching = ((30+31+32)/3) * ones(len(S2branching))

S2massive = 100*array([2.5/5, 1/3.5, 2/5, 1/3, 2/2.5, 1.5/2])
T2massive = ((30+31+32)/3) * ones(len(S2massive))

S2Location = "CAR" #"Florida Keys"

sub2.plot(T2branching-T2_initial, 100-S2branching, "o", markersize = m0, alpha = a0, color = colorbranching)

sub2.plot(T2massive-T2_initial, 100-S2massive, "o", markersize = m0, alpha = a0, color = colormassive)

print(1, T2branching-T2_initial)
print(2,T2massive-T2_initial)

if S2Location not in label_unique and S2Location in Loc_target:   
    sub3.plot(T2-T2_initial, 100-S2, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S2Location)], label = S2Location)
    label_unique.append(S2Location)
elif S2Location in Loc_target:
    sub3.plot(T2-T2_initial, 100-S2, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S2Location)])


#[3] # I remove this because there is not indication of initial temperaturee, nor of MMMax
S3 = array([100 - 86, 100 - 57])
DT3 = ((0.5+1)/2)*ones(len(S3))
#plt.plot(DT3, S3, "o", markersize = m0, label = "[3]")

S3Plate = array([100 - 57])
DT3Plate = ((0.5+1)/2) * ones(len(S3Plate))

S3massive = array([100 - 86])
DT3massive = ((0.5+1)/2) * ones(len(S3massive))

S3Location = "CAR" #"Caribbean"

#sub2.plot(DT3Plate, S3Plate, "o", markersize = m0, color = "grey", label = "Plate-like")
#sub5.plot(DT3massive, S3massive, "o", markersize = m0, color = colormassive)
#sub3.plot(DT3, S3, "o", markersize = m0, label = S3Location)  # I do not include because there is not indication of initial temperature or MMMax
#sub6.plot(DT3, S3, "o", markersize = m0, label = S3Location)  # I do not include because there is not indication of initial temperature or MMMax


#[4]
# temperature not reported

#[5]
T5 = array([30, 32, 30, 32])
T5_initial = 27
T5_MMMax = 30
S5_init = 100*array([0.46/0.5, 0.2/0.5, 1.5/1.75, 0.4/1.75])
S5_MMMax = 100*array([0.46/0.46, 0.2/0.46, 1.5/1.5, 0.4/1.5])

T5branching = array([30, 32, 30, 32])
S5branching_init = 100*array([0.46/0.5, 0.2/0.5, 1.5/1.75, 0.4/1.75])
S5branching_MMMax = 100*array([0.46/0.46, 0.2/0.46, 1.5/1.5, 0.4/1.5])

S5Location = "GBR"

sub2.plot(T5branching - T5_initial, 100-S5branching_init, "o", markersize = m0, alpha = a0, color = colorbranching, label = "Branching")

if S5Location not in label_unique and S5Location in Loc_target:
    sub3.plot(T5-T5_initial, 100-S5_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S5Location)], label = S5Location)
    label_unique.append(S5Location)
elif S5Location in Loc_target:
    sub3.plot(T5-T5_initial, 100-S5_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S5Location)])
    
#[6] 
# Compared to before experiment, since temperature is mean summer all should be MMMax
DT6a_MMMax = array([30-28.8, 30-27.5, 31-28.8, 31-27.5, 32-29.65])
S6a = 100*array([0.75/1.2, 0.9/1.4, 0.40/1.2, 0.2/1.4, 0.5/1.8])

DT6branching_MMMax = array([30-28.8, 30-27.5, 31-28.8, 31-27.5, 32-29.65])
S6branching = 100*array([0.75/1.2, 0.9/1.4, 0.40/1.2, 0.2/1.4, 0.5/1.8])

# Compared to control 
T6b = array([30, 30, 31, 31, 32]) 
T6b_initial = 27.5
S6b_init = 100*array([0.75/0.9, 0.9/1.5, 0.4/0.9, 0.2/1.5, 0.5/1.5])

T6branching2 = array([30, 30, 31, 31, 32])
S6branching2 = 100*array([0.75/0.9, 0.9/1.5, 0.4/0.9, 0.2/1.5, 0.5/1.5])

S6Location = "GBR"
sub2.plot(T6branching2-T6b_initial, 100-S6branching2, "o", markersize = m0, alpha = a0, color = colorbranching)

if S6Location not in label_unique and S6Location in Loc_target:
    sub3.plot(T6b-T6b_initial, 100-S6b_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S6Location)], label = S6Location)
    label_unique.append(S6Location)
elif S6Location in Loc_target:
    sub3.plot(T6b-T6b_initial, 100-S6b_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S6Location)])

#[7]
T7 = array([(33.03 + 33.81)/2])
T7_initial = (28.51 + 29.77)/2
T7_MMMax = 29.77
S7 = array([75])

S7massive = array([75])
T7massive = array([(33.03 + 33.81)/2])

S7Location = "SEA" #"Andaman Sea, Thailand"
sub2.plot(T7massive - T7_initial, 100-S7massive, "o", markersize = m0, alpha = a0, color = colormassive, label = "Massive")

if S7Location not in label_unique and S7Location in Loc_target:
    sub3.plot(T7 - T7_initial, 100-S7, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S7Location)], label = S7Location)
    label_unique.append(S7Location)
elif S7Location in Loc_target:
    sub3.plot(T7 - T7_initial, 100-S7, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S7Location)])

#[8] # This is directly compared with MMMax
S8 = 100 - array([(66+65)/2 , (44+48)/2])
T8 = ((32+34)/2)*ones(len(S8))
T8_MMMax = 30

S8branching = 100 - array([(66+65)/2 , (44+48)/2])
T8branching = ((32+34)/2)*ones(len(S8branching))

S8Location = "GBR"


#[9] 
T9 = array([28, 30])
T9_initial = 24
T9_MMMax = 28
S9_init = 100*array([0.7/1.2, 0.3/1.2])
S9_MMMax = 100*array([0.7/0.7, 0.3/0.7])

#T9massive = array([28, 30]) 
#T9massive_init = 24
#T9massive_MMMax = 28
#S9massive_init = 100*array([0.7/1.2, 0.3/1.2])
#S9massive_MMMax = 100*array([0.7/0.7, 0.3/0.7])

S9Location = "Other"#"Okinawa, Japan"

#sub2.plot(T9massive-T9_initial, 100-S9massive_init, "o", markersize = m0, alpha = a0, color = colormassive)
if S9Location not in label_unique and S9Location in Loc_target:
    sub3.plot(T9 - T9_initial, 100-S9_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S9Location)], label = S9Location)
    label_unique.append(S9Location)
elif S9Location in Loc_target:
    sub3.plot(T9 - T9_initial, 100-S9_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S9Location)])
    
#[10]
S10 = 100*array([0.7/2, 0.7/3, 1.5/3, 1.5/3])
T10 = (30.4)*ones(len(S10))
T10_initial = 26.8
T10_MMMax = 30

S10branching = 100*array([0.7/2, 0.7/3, 1.5/3, 1.5/3])
T10branching = (30.4)*ones(len(S10branching))

S10Location = "Other"#"Mauritus"

#sub2.plot(T10branching - T10_initial, S10branching, "o", markersize = m0, alpha = a0, color = colorbranching)
if S10Location not in label_unique and S10Location in Loc_target:
    sub3.plot(T10 - T10_initial, 100-S10, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S10Location)], label = S10Location)
    label_unique.append(S10Location)
elif S10Location in Loc_target:
    sub3.plot(T10 - T10_initial, 100-S10, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S10Location)])   

#[11]
S11 = 100*array([0.3/1, 0.15/1, 0.4/1, 0.3/5, 0.15/5, 0.4/5])
T11 = (29.75)*ones(len(S11))
T11_initial = 27.5
T11_MMMax = 28.5

S11branching = 100*array([0.3/1, 0.15/1, 0.4/1, 0.3/5, 0.15/5, 0.4/5])
T11branching = (29.75)*ones(len(S11))

S11Location = "Other"#"French Polynesia"
#sub2.plot(T11branching - T11_initial, S11branching, "o", markersize = m0, alpha = a0, color = colorbranching)

if S11Location not in label_unique and S11Location in Loc_target: 
    sub3.plot(T11 - T11_initial, 100-S11, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S11Location)], label = S11Location)
    label_unique.append(S11Location)
elif S11Location in Loc_target:
    sub3.plot(T11 - T11_initial, 100-S11, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S11Location)])

#[12] # directly compare with MMMax because there was no initial temperature to compare with
"""
S12 = 100*array([0.5/1.1, 1/2.3, 0.35/1.2])
T12 = (31.5)*ones(len(S12))
T12_MMMax = 30.5

S12branching = 100*array([0.5/1.1])
T12branching = (31.5)*ones(len(S12branching))

S12massive = 100*array([1/2.3, 0.35/1.2])
T12massive = (31.5)*ones(len(S12massive))

S12Location = "CAR"

sub2.plot(T12branching - T12_MMMax, 100-S12branching, "o", markersize = m0, alpha = a0, color = colorbranching)
sub2.plot(T12massive - T12_MMMax, 10-S12massive, "o", markersize = m0, alpha = a0, color = colormassive)
print(3, T12branching-T12_MMMax)
print(4,T12massive-T12_MMMax)

if S12Location not in label_unique and S12Location in Loc_target:
    sub3.plot(T12 - T12_MMMax, 100-S12, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S12Location)], label = S12Location)
    label_unique.append(S12Location)
elif S12Location in Loc_target:
    sub3.plot(T12 - T12_MMMax, 100-S12, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S12Location)])
"""
#[13]
S13_init = 100*array([0.5/1.6, 0.1/1.6, 0.5/0.95, 0.1/0.95, 0.5/1.1, 0.1/1.1])
T13 = array([32, 34, 32, 34, 32, 34])
T13_initial = 26

T13_MMMax = 30
S13_MMMax = 100*array([0.5/2, 0.1/2, 0.5/0.75, 0.1/0.75, 0.5/0.8, 0.1/0.8])

T13massive = array([32, 34])
S13massive_init = 100*array([0.5/1.6, 0.1/1.6])
S13massive_MMMax = 100*array([0.5/2, 0.1/2]) # compared with MMMax, not considered

T13otherShape = array([32, 34, 32, 34])
S13otherShape_init = 100*array([0.5/0.95, 0.1/0.95, 0.5/1.1, 0.1/1.1])
S13otherShape_MMMax = 100*array([0.5/0.75, 0.1/0.75, 0.5/0.8, 0.1/0.8]) # compared with MMMax, not considered

S13Location = "CAR"
sub2.plot(T13massive - T13_initial, 100-S13massive_init, "o", markersize = m0, alpha = a0, color = colormassive)
sub2.plot(T13otherShape - T13_initial, 100-S13otherShape_init, "o", markersize = m0, alpha = a0, color = colorotherShape, label = "Others")

print(5, T13massive-T13_initial)
print(6,T13otherShape-T13_initial)

if S13Location not in label_unique and S13Location in Loc_target:
    sub3.plot(T13 - T13_initial, 100-S13_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S13Location)], label = S13Location)
    label_unique.append(S13Location)
elif S13Location in Loc_target:
    sub3.plot(T13 - T13_initial, 100-S13_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S13Location)])

#[14]
S14_init = 100*array([4.31/4.71, 1.77/4.71, 0.135/4.71, 7.77/10.02, 4.68/10.02])
T14 = array([28.44, 29.61, 31.68, 27.89, 30.37])
T14_initial = array([27.87, 27.87, 27.87, 26.21, 26.21])

S14_MMMax = 100*array([4.31/4.31, 1.77/4.31, 0.135/4.31, 7.77/7.77, 4.68/7.77])
T14_MMMax = array([28.44, 28.44, 28.44, 27.89, 27.89])
DTMMMax = T14 - T14_MMMax

S14branching_init = S14_init
T14branching = T14
S14branching_MMMax = S14_MMMax
T14branching = T14

S14Location = "CAR"
sub2.plot(T14branching - T14_initial, 100-S14branching_init, "o", markersize = m0, alpha = a0, color = colorbranching)
DTmmmax = T14branching - T14_MMMax

print(7, T14branching - T14_initial)

if S14Location not in label_unique and S14Location in Loc_target:
    sub3.plot(T14 - T14_initial, 100-S14_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S14Location)], label = S14Location)
    label_unique.append(S14Location)
elif S14Location in Loc_target:
    sub3.plot(T14 - T14_initial, 100-S14_init, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S14Location)])

#[15]
S15 = 100*array([4./8, 5.5/16, 6./14, 7./15, 6.5/11])
T15 = (30.5)*ones(len(S15))
T15_initial = (27.+29.)/2
T15_MMMax = 29.

S15branching = 100*array([4./8, 5.5/16]) 
T15branching = 30.5*ones(len(S15branching))

S15massive = 100*array([6./14, 7./15, 6.5/11])
T15massive = 30.5*ones(len(S15massive))

S15Location = "CAR" 
sub2.plot(T15branching - T15_initial, 100-S15branching, "o", markersize = m0, alpha = a0, color = colorbranching)
sub2.plot(T15massive - T15_initial, 100-S15massive, "o", markersize = m0, alpha = a0, color = colormassive)

print(8, T15branching - T15_initial)
print(8, T15massive - T15_initial)

if S15Location not in label_unique and S15Location in Loc_target:
    sub3.plot(T15 - T15_initial, 100-S15, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S15Location)], label = S15Location)
    label_unique.append(S15Location)
elif S15Location in Loc_target:
    sub3.plot(T15 - T15_initial, 100-S15, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S15Location)])

#[16]
S16 = 100*array([0.01/1, 0.01/10, 0.01/10, 0.5/10, 0.8/10, 0.7/5])
T16 = array([30.5, 30.5, 30.5, 30, 30, 30])
T16_initial = 28.

T16_MMMax = 29

S16branching = 100*array([0.01/10, 0.01/10, 0.5/10, 0.8/10])
T16branching = array([30.5, 30.5, 30, 30])

S16massive = 100*array([0.01/1, 0.7/5])
T16massive = array([30.5, 30])

sub2.plot(T16branching - T16_initial, 100-S16branching, "o", markersize = m0, alpha = a0, color = colorbranching)
sub2.plot(T16massive - T16_initial, 100-S16massive, "o", markersize = m0, alpha = a0, color = colormassive)


S16Location = "CAR"

print(10, T16branching - T16_initial)
print(11, T15massive - T15_initial)
    

if S16Location not in label_unique and S16Location in Loc_target:
    sub3.plot(T16 - T16_initial, 100-S16, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S16Location)], label = S16Location)
    label_unique.append(S16Location)
elif S16Location in Loc_target:
    sub3.plot(T16 - T16_initial, 100-S16, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S16Location)])

#[17] 

S17 = 100*array([3./10, 6./11, 1.5/9, 4./18, 1./18, 8./12])
T17 = 30.25*ones(len(S17))
T17_initial = 28.5

T17_MMMax = 29.5

S17massive = S17
T17massive = T17
sub2.plot(T17massive - T17_initial, 100-S17massive, "o", markersize = m0, alpha = a0, color = colormassive)

S17Location = "SEA"

if S17Location not in label_unique and S17Location in Loc_target:
    sub3.plot(T17 - T17_initial, 100-S17, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S17Location)], label = S17Location)
    label_unique.append(S17Location)
elif S17Location in Loc_target:
    sub3.plot(T17 - T17_initial, 100-S17, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S17Location)])    

#[18] No MMMax

S18 = 100*array([5.52/7.4])
T18 = 28.81
T18_initial = 28.18

S18massive = S18
T18massive = T18

sub2.plot(T18massive - T18_initial, 100-S18massive, "o", markersize = m0, alpha = a0, color = colormassive)

S18Location = "SEA"

if S18Location not in label_unique and S18Location in Loc_target:
    sub3.plot(T18 - T18_initial, 100-S18, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S18Location)], label = S18Location)
    label_unique.append(S18Location)
elif S17Location in Loc_target:
    sub3.plot(T18 - T18_initial, 100-S18, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(S18Location)])


sub2.legend(ncol = 1, fontsize = fsize)

sub3.legend(ncol = 1, fontsize = fsize)

plt.savefig("Figures/EPS/FigS1a.eps", dpi= 600, bbox_inches = 'tight') 
plt.savefig("Figures/PDF/FigS1a.pdf", dpi= 600, bbox_inches = 'tight')
    

                                                                        
""" Separate by region """
fsize = 18*3.25#20
fsize2 = 16*3.25#13
plt.rcParams["font.family"] = "arial"   

font0 = FontProperties()
font = font0.copy()
weight = "bold"
font.set_weight(weight) 

fig = plt.figure(figsize = (50, 14))
sub2 = fig.add_subplot(1, 3, 1)
plt.title("GBR", fontsize = fsize, fontproperties = font)
plt.ylabel("Reduction in symbiont biomass \n $\Delta S$ (%)", fontsize = fsize)
plt.xlabel("$\Delta T$"+u" (\N{DEGREE SIGN}C)", fontsize = fsize)
plt.text(-1.5, 101, "c", fontproperties = font, fontsize = fsize2)

plt.xticks((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), fontsize = fsize2)
#plt.text(-1.2, 100, "A", fontsize = fsize, fontproperties = font)
plt.xlim(0, 9)

sub2.tick_params(axis = "y", direction = "in", labelsize = fsize2)
#sub2.set_yticks((0, 25, 50, 75, 100)) # GBR and SEA
ylabs = arange(0, 110, 10)
sub2.set_yticks(ylabs) # GBR and SEA
sub2.set_yticklabels(list(["%d"%s for s in ylabs]))
sub2.set_ylim((0, 100))
part2 = sub2.twinx()
part2.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
part2.set_yticks(ylabs)
part2.set_yticklabels([" " for d in ylabs])
part2.set_ylim((0, 100))


sub3 = fig.add_subplot(1, 3, 2)
plt.title("SEA", fontsize = fsize, fontproperties = font)
plt.xlabel("$\Delta T$"+u" (\N{DEGREE SIGN}C)", fontsize = fsize)
plt.xticks((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), fontsize = fsize2)
plt.xlim(0, 9)
plt.text(-0.5, 101, "d", fontproperties = font, fontsize = fsize2)

sub3.tick_params(axis = "y", direction = "in", labelsize = fsize2)
sub3.set_yticks(ylabs) # GBR and SEA
sub3.set_yticklabels(list([" "%s for s in ylabs]))
sub3.set_ylim((0, 100))
part3 = sub3.twinx()
part3.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
part3.set_yticks(ylabs)
part3.set_yticklabels([" " for d in ylabs])
part3.set_ylim((0, 100))

sub4 = fig.add_subplot(1, 3, 3)
plt.title("CAR", fontsize = fsize, fontproperties = font)
plt.xlabel("$\Delta T$"+u" (\N{DEGREE SIGN}C)", fontsize = fsize)
plt.text(-0.5, 101, "e", fontproperties = font, fontsize = fsize2)

plt.xticks((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), fontsize = fsize2)
plt.xlim(0, 9)
sub4.tick_params(axis = "y", direction = "in", labelsize = fsize2)
sub4.set_yticks(ylabs) # GBR and SEA
sub4.set_yticklabels(list([" "%s for s in ylabs]))
sub4.set_ylim((0, 100))
part4 = sub4.twinx()
part4.tick_params(axis = "y", direction = "in", labelsize = fsize2) 
part4.set_yticks(ylabs)
part4.set_yticklabels([" " for d in ylabs])
part4.set_ylim((0, 100))


a0 = 1
m0 = 20
def plot_func(T, T_init, S, Loc):
    if Loc == "GBR":
        sub2.plot(T-T_init, 100-S, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(Loc)], label = Loc)
    elif Loc == "SEA":
        sub3.plot(T-T_init, 100-S, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(Loc)], label = Loc)
    elif Loc == "CAR":
        sub4.plot(T-T_init, 100-S, "o", markersize = m0, alpha = a0, color = col[Loc_target.index(Loc)], label = Loc)
    

#color for references
numLoc = 16
cmx = mpl.cm
cmapC = mpl.cm.tab20c#Vega20c
cNorm = mpl.colors.Normalize(vmin=1, vmax=numLoc)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmapC)

colorbranching = "blue"
colormassive = "orange"
colorotherShape = "grey"

Loc_target = ["GBR", "SEA", "CAR"]
#col = [(0.827, 0.553, 0.686), (0.416, 0.071, 0.239), (0.69, 0.345, 0.514), "grey"]
col = ["green",(0.44, 0.07, 0.20), (0.8, 0.4, 0.8), (0.9, 0.9, 0.9)]  
label_unique = [] # for regions

#[1] separate between those with visual signs and those with no visual signs of bleaching
S1 = array([(10+11)/2, (28+16)/2, (52+47)/2, (55+69)/2])
T1 = (30.8)*ones(len(S1))
T1_initial = 28.6  
T1_MMMax = 29.8

S1branching = array([(10+11)/2, (52+47)/2])
T1branching = (30.8)*ones(len(S1branching))

S1massive = array([(28+16)/2, (55+69)/2])
T1massive = (30.8)*ones(len(S1massive))

S1Location = "Other"#"South China Sea"

plot_func(T1, T1_initial, S1, S1Location)
    
#[2]
S2 = 100*array([1.5/3.5, 0.9/1.5, 2.5/5, 1/3.5, 2/5, 1/3, 1/1.5, 2/2.5, 1.5/2])
T2 = ((30+31+32)/3) * ones(len(S2))
T2_initial = (26 + 27)/2
T2_MMMax = 29

S2branching = 100*array([1.5/3.5, 0.9/1.5, 1/1.5])
T2branching = ((30+31+32)/3) * ones(len(S2branching))

S2massive = 100*array([2.5/5, 1/3.5, 2/5, 1/3, 2/2.5, 1.5/2])
T2massive = ((30+31+32)/3) * ones(len(S2massive))

S2Location = "CAR" #"Florida Keys"

plot_func(T2branching, T2_initial, S2branching, S2Location)
plot_func(T2massive, T2_initial, S2massive, S2Location)


#[3] # I remove this because there is not indication of initial temperaturee, nor of MMMax
S3 = array([100 - 86, 100 - 57])
DT3 = ((0.5+1)/2)*ones(len(S3))
#plt.plot(DT3, S3, "o", markersize = m0, label = "[3]")

S3Plate = array([100 - 57])
DT3Plate = ((0.5+1)/2) * ones(len(S3Plate))

S3massive = array([100 - 86])
DT3massive = ((0.5+1)/2) * ones(len(S3massive))

S3Location = "CAR" #"Caribbean"


#[4]
# temperature not reported

#[5]
T5 = array([30, 32, 30, 32])
T5_initial = 27
T5_MMMax = 30
S5_init = 100*array([0.46/0.5, 0.2/0.5, 1.5/1.75, 0.4/1.75])
S5_MMMax = 100*array([0.46/0.46, 0.2/0.46, 1.5/1.5, 0.4/1.5])

T5branching = array([30, 32, 30, 32])
S5branching_init = 100*array([0.46/0.5, 0.2/0.5, 1.5/1.75, 0.4/1.75])
S5branching_MMMax = 100*array([0.46/0.46, 0.2/0.46, 1.5/1.5, 0.4/1.5])

S5Location = "GBR"

plot_func(T5branching, T5_initial, S5branching_init, S5Location)

    
#[6] 
# Compared to before experiment, since temperature is mean summer all should be MMMax
DT6a_MMMax = array([30-28.8, 30-27.5, 31-28.8, 31-27.5, 32-29.65])
S6a = 100*array([0.75/1.2, 0.9/1.4, 0.40/1.2, 0.2/1.4, 0.5/1.8])

DT6branching_MMMax = array([30-28.8, 30-27.5, 31-28.8, 31-27.5, 32-29.65])
S6branching = 100*array([0.75/1.2, 0.9/1.4, 0.40/1.2, 0.2/1.4, 0.5/1.8])

# Compared to control 
T6b = array([30, 30, 31, 31, 32]) 
T6b_initial = 27.5
S6b_init = 100*array([0.75/0.9, 0.9/1.5, 0.4/0.9, 0.2/1.5, 0.5/1.5])

T6branching2 = array([30, 30, 31, 31, 32])
S6branching2 = 100*array([0.75/0.9, 0.9/1.5, 0.4/0.9, 0.2/1.5, 0.5/1.5])

S6Location = "GBR"

plot_func(T6branching2, T6b_initial, S6branching2, S6Location)

#[7]
T7 = array([(33.03 + 33.81)/2])
T7_initial = (28.51 + 29.77)/2
T7_MMMax = 29.77
S7 = array([75])

S7massive = array([75])
T7massive = array([(33.03 + 33.81)/2])

S7Location = "SEA" #"Andaman Sea, Thailand"

plot_func(T7massive, T7_initial, S7massive, S7Location)


#[8] # This is directly compared with MMMax
S8 = 100 - array([(66+65)/2 , (44+48)/2])
T8 = ((32+34)/2)*ones(len(S8))
T8_MMMax = 30

S8branching = 100 - array([(66+65)/2 , (44+48)/2])
T8branching = ((32+34)/2)*ones(len(S8branching))

S8Location = "GBR"


#[9] 

T9massive = array([28, 30]) 
T9massive_init = 24
T9massive_MMMax = 28
S9massive_init = 100*array([0.7/1.2, 0.3/1.2])
S9massive_MMMax = 100*array([0.7/0.7, 0.3/0.7])

S9Location = "Other"#"Okinawa, Japan"

plot_func(T9massive, T9massive_init, S9massive_init, S9Location)


#[10]
T10_initial = 26.8
T10_MMMax = 30

S10branching = 100*array([0.7/2, 0.7/3, 1.5/3, 1.5/3])
T10branching = (30.4)*ones(len(S10branching))

S10Location = "Other"#"Mauritus"

plot_func(T10branching, T10_initial, S10branching, S10Location)
 

#[11]
S11 = 100*array([0.3/1, 0.15/1, 0.4/1, 0.3/5, 0.15/5, 0.4/5])
T11 = (29.75)*ones(len(S11))
T11_initial = 27.5
T11_MMMax = 28.5

S11branching = 100*array([0.3/1, 0.15/1, 0.4/1, 0.3/5, 0.15/5, 0.4/5])
T11branching = (29.75)*ones(len(S11))

S11Location = "Other"#"French Polynesia"
#sub2.plot(T11branching - T11_initial, S11branching, "o", markersize = m0, alpha = a0, color = colorbranching)
plot_func(T11branching, T11_initial, 100-S11branching, S11Location)

#[12] # directly compare with MMMax because there was no initial temperature to compare with
"""
S12 = 100*array([0.5/1.1, 1/2.3, 0.35/1.2])
T12 = (31.5)*ones(len(S12))
T12_MMMax = 30.5

S12branching = 100*array([0.5/1.1])
T12branching = (31.5)*ones(len(S12branching))

S12massive = 100*array([1/2.3, 0.35/1.2])
T12massive = (31.5)*ones(len(S12massive))

S12Location = "CAR"

plot_func(T12branching, T12_MMMax, 100-S12branching, S12Location)
plot_func(T12massive, T12_MMMax, 100-S12massive, S12Location)
"""

#[13]
S13_init = 100*array([0.5/1.6, 0.1/1.6, 0.5/0.95, 0.1/0.95, 0.5/1.1, 0.1/1.1])
T13 = array([32, 34, 32, 34, 32, 34])
T13_initial = 26

T13_MMMax = 30
S13_MMMax = 100*array([0.5/2, 0.1/2, 0.5/0.75, 0.1/0.75, 0.5/0.8, 0.1/0.8])

T13massive = array([32, 34])
S13massive_init = 100*array([0.5/1.6, 0.1/1.6])
S13massive_MMMax = 100*array([0.5/2, 0.1/2]) # compared with MMMax, not considered

T13otherShape = array([32, 34, 32, 34])
S13otherShape_init = 100*array([0.5/0.95, 0.1/0.95, 0.5/1.1, 0.1/1.1])
S13otherShape_MMMax = 100*array([0.5/0.75, 0.1/0.75, 0.5/0.8, 0.1/0.8]) # compared with MMMax, not considered

S13Location = "CAR"

plot_func(T13massive, T13_initial,S13massive_init, S13Location)
plot_func(T13otherShape, T13_initial,S13otherShape_init, S13Location)


#[14]
S14_init = 100*array([4.31/4.71, 1.77/4.71, 0.135/4.71, 7.77/10.02, 4.68/10.02])
T14 = array([28.44, 29.61, 31.68, 27.89, 30.37])
T14_initial = array([27.87, 27.87, 27.87, 26.21, 26.21])

S14_MMMax = 100*array([4.31/4.31, 1.77/4.31, 0.135/4.31, 7.77/7.77, 4.68/7.77])
T14_MMMax = array([28.44, 28.44, 28.44, 27.89, 27.89])
DTMMMax = T14 - T14_MMMax

S14branching_init = S14_init
T14branching = T14
S14branching_MMMax = S14_MMMax
T14branching = T14

S14Location = "CAR"
DTmmmax = T14branching - T14_MMMax
plot_func(T14branching, T14_initial, S14branching_init, S14Location)

#[15]
S15 = 100*array([4./8, 5.5/16, 6./14, 7./15, 6.5/11])
T15 = (30.5)*ones(len(S15))
T15_initial = (27.+29.)/2
T15_MMMax = 29.

S15branching = 100*array([4./8, 5.5/16]) 
T15branching = 30.5*ones(len(S15branching))

S15massive = 100*array([6./14, 7./15, 6.5/11])
T15massive = 30.5*ones(len(S15massive))

S15Location = "CAR" 

plot_func(T15branching, T15_initial, S15branching, S15Location)
plot_func(T15massive, T15_initial, S15massive, S15Location)


#[16]
S16 = 100*array([0.01/1, 0.01/10, 0.01/10, 0.5/10, 0.8/10, 0.7/5])
T16 = array([30.5, 30.5, 30.5, 30, 30, 30])
T16_initial = 28.

T16_MMMax = 29

S16branching = 100*array([0.01/10, 0.01/10, 0.5/10, 0.8/10])
T16branching = array([30.5, 30.5, 30, 30])

S16massive = 100*array([0.01/1, 0.7/5])
T16massive = array([30.5, 30])

S16Location = "CAR"
plot_func(T16branching, T16_initial, S16branching, S16Location)
plot_func(T16massive, T16_initial, S16massive, S16Location)

#[17] 

S17 = 100*array([3./10, 6./11, 1.5/9, 4./18, 1./18, 8./12])
T17 = 30.25*ones(len(S17))
T17_initial = 28.5

T17_MMMax = 29.5

S17massive = S17
T17massive = T17

S17Location = "SEA"

plot_func(T17massive, T17_initial, S17massive, S17Location)


#[18] No MMMax

S18 = 100*array([5.52/7.4])
T18 = 28.81
T18_initial = 28.18

S18massive = S18
T18massive = T18

S18Location = "SEA"

plot_func(T18massive, T18_initial, S18massive, S18Location)


plt.savefig("Figures/EPS/FigS1b.eps", dpi= 600, bbox_inches = 'tight') 
plt.savefig("Figures/PDF/FigS1b.pdf", dpi= 600, bbox_inches = 'tight')

                                                                                                                                                                                                                                                                                 

