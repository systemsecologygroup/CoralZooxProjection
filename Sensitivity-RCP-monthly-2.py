# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, arange, isnan
from scipy.integrate import ode
from matplotlib import pyplot as plt
import matplotlib
plt.switch_backend('agg')

import pdb

"""Generates data for Figure 5 & 6"""

#### Model sensitivity for simulations with bleaching, time in month #####
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

gammaH = 1e6 # free param 
fsize = 18

r = (12)*1e3
error = 1e-2 # prevent division by zeros these scales are alright because we are dealing with very high numbers
error2 = 1e-2 # prevent division by zeros
error3 = 1e-2 # prevent division by zeros


def Gradient(coral, u, symbiont, Gx1, beta, alpha, K_symb, K_C_Reg):
    symbiontH = gammaH*coral
    kappa = symbiont/(symbiontH + symbiont + error3) 
    cost_gam = symbiont/(K_symb + error)
    Grad = Gx1*beta*kappa*(1 - coral/K_C_Reg)*exp(-beta*u) - r*alpha*cost_gam*exp(r*u)
    return Grad

# Future of coral biomass according to IPCC emission senarios B1 
AddTime = 2000*12 # for spine up, longer time to reach a steady state
def SystemForcing(t, y, T0, rho, skew, N, NormCor, TempNS, K_C_Reg):
    dSystem = zeros(len(y))
    coral = y[0]
    u = y[1]
    E = 1 - exp(-beta*u)
    #E = max(0, E) # make sure E is not above 1, due to fast decrease of u making u<0
    symbiont = y[2]
    # Introducing temperature dependence
    #print t
    if t<=NumTime+AddTime:
        Temperature1 = TempNS[int(t)] # this should rather be some function of t and not index but it works here because simulation step = 1
    else:
        Temperature1 = TempNS[NumTime+AddTime] # there are index out of range sometimes, maybe because of the integrator, so I just make sure not to get those errors
    Tcenter = (Temperature1-T0)/rho
    Gx1Forcing = G_C*norm.pdf(Tcenter)*norm.cdf(Tcenter*skew)/max(NormCor) 
    # Computing the derivative
    K_symb = Ksmax*coral 
    symbiontH = gammaH*coral
    kappa = symbiont/(symbiontH + symbiont + error3)
    cost_gam = symbiont/(K_symb + error)
    Benefit = Gx1Forcing*kappa*E*(1-coral/K_C_Reg)
    Cost =alpha*exp(r*u)*(cost_gam) + M_C
    Fitness = (Benefit - Cost)
    dSystem[0] = Fitness*coral
    dSystem[1] = N*Gradient(coral, u, symbiont, Gx1Forcing, beta, alpha, K_symb, K_C_Reg)
    G_S =  a*exp(b*Temperature1) # symbiont temperature associated growth
    dSystem[2] =G_S*(1-symbiont/(K_symb + error2))*symbiont
    #print K_symb, coral, symbiont, dSystem
    return dSystem

# Parameter resp. for the different region GBR, SEA , CAR
T0_list = array([26.78, 28.13, 27.10]) 
skew_list = array([0.0002, 3.82, 1.06])
rho_list = array([1.0, 0.81, 0.89])

color_list = array([(0.647*col, 0.333*col, 0.075*col) for col in xrange(2,5+2)])/5

scale = 1e-11 

rawNum_GBR_orig = arange(0.01, 0.1, 0.00005)
rawNum_SEA_orig = arange(0.01, 0.1, 0.00005)
rawNum_CAR_orig = arange(0.01, 0.1, 0.00005) 



# This is just to reduce the number of simulation, I did not have time to redo all of everything with the newly estimated values of N, 
# These are not the index of estimated N reported in manuscript, it was from earlier code versions, the choice does not matter
GBR_N_index = 100 # index in rawNum_GBR  0.0145
SEA_N_index = 275 # index in rawNum_SEA  0.02375
CAR_N_index = 250 # index in rawNum_CAR  0.0225
"""
# Only for sensitivity simulations to reduce the number of simulations (USED for supplementary figure in manuscript)
rawNum_GBR = concatenate((arange(min(rawNum_GBR_orig), rawNum_GBR_orig[GBR_N_index], 0.0005), array([0.75*rawNum_GBR_orig[GBR_N_index]]), array([rawNum_GBR_orig[GBR_N_index]]), array([1.25*rawNum_GBR_orig[GBR_N_index]]), arange(rawNum_GBR_orig[GBR_N_index]+0.0005, max(rawNum_GBR_orig), 0.0005)))
rawNum_SEA = concatenate((arange(min(rawNum_SEA_orig), rawNum_SEA_orig[SEA_N_index], 0.0005), array([0.75*rawNum_SEA_orig[SEA_N_index]]), array([rawNum_SEA_orig[SEA_N_index]]), array([1.25*rawNum_SEA_orig[SEA_N_index]]), arange(rawNum_SEA_orig[SEA_N_index]+0.0005, max(rawNum_SEA_orig), 0.0005)))
rawNum_CAR = concatenate((arange(min(rawNum_CAR_orig), rawNum_CAR_orig[CAR_N_index], 0.0005), array([0.75*rawNum_CAR_orig[CAR_N_index]]), array([rawNum_CAR_orig[CAR_N_index]]), array([1.25*rawNum_CAR_orig[CAR_N_index]]), arange(rawNum_CAR_orig[CAR_N_index]+0.0005, max(rawNum_CAR_orig), 0.0005))) 
"""
# Only for sensitivity simulations to reduce the number of simulations (Just for a DEMO here)
rawNum_GBR = concatenate((arange(min(rawNum_GBR_orig), rawNum_GBR_orig[GBR_N_index], 0.05), array([rawNum_GBR_orig[GBR_N_index]]), arange(rawNum_GBR_orig[GBR_N_index]+0.05, max(rawNum_GBR_orig), 0.05)))
rawNum_SEA = concatenate((arange(min(rawNum_SEA_orig), rawNum_SEA_orig[SEA_N_index], 0.05), array([rawNum_SEA_orig[SEA_N_index]]), arange(rawNum_SEA_orig[SEA_N_index]+0.05, max(rawNum_SEA_orig), 0.05)))
rawNum_CAR = concatenate((arange(min(rawNum_CAR_orig), rawNum_CAR_orig[CAR_N_index], 0.05), array([rawNum_CAR_orig[CAR_N_index]]), arange(rawNum_CAR_orig[CAR_N_index]+0.05, max(rawNum_CAR_orig), 0.05))) 

# Parameters to study sensitivity to
ParamsPlus= ["Gmax-Plus", "a-Plus", "b-Plus", "KC-Plus", "Ksmax-Plus", "MC-Plus", "alpha-Plus", "r-Plus", "beta-Plus", "GammaH-Plus"]
ParamsMinus= ["Gmax-Minus", "a-Minus", "b-Minus", "KC-Minus", "Ksmax-Minus", "MC-Minus", "alpha-Minus", "r-Minus", "beta-Minus", "GammaH-Minus"]
Params = ParamsPlus + ParamsMinus


# Fixed control parameters
G_C0 = G_C
a0 = a
b0 = b
alpha0 = alpha
M_C0 = M_C
K_C0 = K_C_List
Ksmax0 = Ksmax
beta0 = beta
gammaH0 = gammaH
r0 = r

TempList = linspace(0, 75, 1500)
RCP_list = ["RCP26", "RCP45", "RCP85"]

Locations = array(["GBR", "SEA", "CAR"])
filename = "Sensitivity-N/Monthly/"
for ind in xrange(len(Params)):
    PlusOrMinus = (-1)*(Params[ind] in ParamsMinus) + (1)*(Params[ind] in ParamsPlus) # only one of the statement is true because ParamsPlus and ParamsMinus have no intersection
    Bool = 0.25*PlusOrMinus
    print ind, PlusOrMinus
    
    Test_if_G_C = (Params[ind] == "Gmax-Plus") + (Params[ind] == "Gmax-Minus")  # will always give 0 or 1
    Test_if_a = (Params[ind] == "a-Plus") + (Params[ind] == "a-Minus")
    Test_if_b = (Params[ind] == "b-Plus") + (Params[ind] == "b-Minus")
    Test_if_alpha = (Params[ind] == "alpha-Plus") + (Params[ind] == "alpha-Minus")
    Test_if_M_C = (Params[ind] == "MC-Plus") + (Params[ind] == "MC-Minus")
    Test_if_K_C = (Params[ind] == "KC-Plus") + (Params[ind] == "KC-Minus")
    Test_if_Ksmax = (Params[ind] == "Ksmax-Plus") + (Params[ind] == "Ksmax-Minus")
    Test_if_beta = (Params[ind] == "beta-Plus") + (Params[ind] == "beta-Minus")
    Test_if_gammaH = (Params[ind] == "GammaH-Plus") + (Params[ind] == "GammaH-Minus")
    Test_if_r = (Params[ind] == "r-Plus") + (Params[ind] == "r-Minus")
    
    G_C = G_C0 + Bool*G_C0*Test_if_G_C  # this is G_max in the model
    a = a0 + Bool*a0*Test_if_a # Symbiont specific growth rate, linear growth rate
    b = b0 + Bool*b0*Test_if_b # Exponential Growth constant of the symbiont
    alpha = alpha0 + Bool*alpha0*Test_if_alpha 
    M_C = M_C0 + Bool*M_C0*Test_if_M_C # coral mortality 
    Ksmax = Ksmax0 + Bool*Ksmax0*Test_if_Ksmax # healthy measure of carying capacity of symbiont per host biomass Gamma
    beta = beta0 + Bool*beta0*Test_if_beta
    gammaH = gammaH0 + Bool*gammaH0*Test_if_gammaH
    r = r0 + Bool*r0*Test_if_r
    
    for RCP in RCP_list:
        rcp = RCP
        for z in xrange(len(Locations)):
            Location_value = Locations[z]

            plt.figure()
            sub1 = plt.subplot(3, 1, 1)
            sub2 = plt.subplot(3, 1, 2)
            sub3 = plt.subplot(3, 1, 3)
            
            #variable time
            file0 = open("Monthly-SST-scenarios/Months-"+Location_value+"-"+rcp+"-MPI"+".dat", "r")
            time = load(file0, allow_pickle = True)   
            file0.close()
            
            #Scenarios
            file1 = open("Monthly-SST-scenarios/SST-"+Location_value+"-"+rcp+"-MPI"+".dat", "r")
            TempNSdata = load(file1, allow_pickle = True)
            file1.close()
            
            TimeMonth = arange(0, len(time), 1.0)
            NumTime = len(TimeMonth) - 1
            
            HOST = []
            TRAIT = []
            SYMB = []
            T0 = T0_list[z]
            skew = skew_list[z]
            rho = rho_list[z]
            K_C_Reg = K_C0[z] + Bool*K_C0[z]*Test_if_K_C # Regional coral carying capacity
            
            TempNS = concatenate((mean(TempNSdata[time<=2000+11/12])*ones(AddTime+12), TempNSdata)) # use mean(TempNSdata[time<=2000]) = mean WOD13 before year 2000 as initial temperature profile for spine up
            TempCORListCenter = (TempList - T0)/rho
            NormCor = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)
            Dist = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)/max(NormCor)
            maxD = max(Dist)
            T_opt = TempList[int(list(Dist).index(maxD))] 
            
            fileCoral = open(filename+"/CORAL-"+RCP+"-"+Locations[z]+"-"+Params[ind]+".dat", "wr")
            fileTrait = open(filename+"/TRAIT-"+RCP+"-"+Locations[z]+"-"+Params[ind]+".dat", "wr")
            fileSymb = open(filename+"/SYMB-"+RCP+"-"+Locations[z]+"-"+Params[ind]+".dat", "wr")
            if z == 0:
                N_List = (scale)*rawNum_GBR
                Initial = array([0.75*K_C_GBR, 0.0000005, 0.001]) 
            elif z == 1:
                N_List = (scale)*rawNum_SEA
                Initial = array([0.75*K_C_SEA, 0.0000005, 0.001]) 
            else:
                N_List = (scale)*rawNum_CAR
                Initial = array([0.75*K_C_CAR, 0.0000005, 0.001]) 
            for i in xrange(len(N_List)):
                N = N_List[i]
                ode15s = ode(SystemForcing)
                ode15s.set_f_params(T0, rho, skew, N, NormCor, TempNS, K_C_Reg)
                ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
                ode15s.set_initial_value(Initial, 0)
                Dynamics = zeros((len(time)+AddTime+12, 3))
                Dynamics[0, 0] = ode15s.y[0]
                Dynamics[0, 1] = ode15s.y[1]
                Dynamics[0, 2] = ode15s.y[2]
                k=1
                while ode15s.successful() and ode15s.t <= len(time)-2+AddTime+12:
                    if Location_value == "GBR":   # bleaching thresholds for the GBR
                        if TempNS[int(ode15s.t)] >= T_opt+2 and TempNS[int(ode15s.t)] < T_opt+4: # bleaching occurs symbionts are expeled and the system is reset to restart again
                            random_reduction = random.uniform(1-0.95, 1-0.10)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1) # simulation step =1
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "reset"
                        elif TempNS[int(ode15s.t)] >= T_opt+4:
                            random_reduction = random.uniform(1-0.95, 1-0.60)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                        else:
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "no reset"
                    elif Location_value == "SEA": # bleaching thresholds for the SEA
                        if TempNS[int(ode15s.t)] >= T_opt+1: # bleaching occurs symbionts are expeled and the system is reset to restart again
                            random_reduction = random.uniform(1-0.95, 1-0.25)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "reset"
                        else:
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "no reset"
                    elif Location_value == "CAR": # bleaching thresholds for the CAR
                        if TempNS[int(ode15s.t)] >= T_opt+1 and TempNS[int(ode15s.t)] < T_opt+3: # bleaching occurs symbionts are expeled and the system is reset to restart again
                            random_reduction = random.uniform(1-0.95, 1-0.15)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "reset"
                        elif TempNS[int(ode15s.t)] >= T_opt+3 and TempNS[int(ode15s.t)] < T_opt+6:
                            random_reduction = random.uniform(1-0.95, 1-0.35)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                        elif TempNS[int(ode15s.t)] >= T_opt+6:
                            random_reduction = random.uniform(1-0.95, 1-0.80)
                            ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                        else:
                            ode15s.integrate(ode15s.t+1)
                            Dynamics[k, 0] = ode15s.y[0]
                            Dynamics[k, 1] = ode15s.y[1]
                            Dynamics[k, 2] = ode15s.y[2]
                            #print rcp, Location_value, i, ode15s.t, "no reset"
                    k+=1
                HOST.append(Dynamics[:, 0])
                TRAIT.append(Dynamics[:, 1])
                SYMB.append(Dynamics[:, 2])
            array(HOST).dump(fileCoral)
            array(TRAIT).dump(fileTrait)
            array(SYMB).dump(fileSymb)
            fileCoral.flush()
            fileTrait.flush()
            fileSymb.flush()
            fileCoral.close()
            fileTrait.close()