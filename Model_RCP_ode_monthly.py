# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, arange, isnan
from scipy.integrate import ode
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import pdb
"""
This code is generates our model dynamics
"""

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

Initial_GBR = array([0.75*K_C_GBR, 0.0000005, 0.001]) 
Initial_SEA = array([0.75*K_C_SEA, 0.0000005, 0.001]) 
Initial_CAR = array([0.75*K_C_CAR, 0.0000005, 0.001]) 

Initial_list = {"GBR":Initial_GBR, "SEA":Initial_SEA, "CAR":Initial_CAR}


Ksmax = 3e6 # healthy measure of carying capacity of symbiont per host biomass Gamma

beta = (12)*1e2

gammaH = 1e6 # minimun symbiont density found on healthy coral colonie
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

# The following function defines the system of ordinary differential equation to be integrated 
AddTime = 2000*12 # for spine up to reach a quasi steady state
def SystemForcing(t, y, T0, rho, skew, N, NormCor, TempNS, K_C_Reg, NumTime):
    dSystem = zeros(len(y))
    coral = y[0]
    u = y[1]
    E = 1 - exp(-beta*u)
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



scale = 1e-11
# Parameter resp. for the different region GBR, SEA , CAR
T0_list = array([26.78, 28.13, 27.10]) 
skew_list = array([0.0002, 3.82, 1.06])
rho_list = array([1.0, 0.81, 0.89])

TempList = linspace(0, 75, 1500)

color_list = array([(0.647*col, 0.333*col, 0.075*col) for col in xrange(2,5+2)])/5

def RUN_SIM(RCP_list, Location_value, N_values, folder, Init = Initial_list):
    if Location_value == "GBR":
        z = 0
    elif Location_value == "SEA":
        z = 1
    elif Location_value == "CAR":
        z = 2
    
    HostDic = {}
    TraitDic = {}
    SymbDic = {}
    
    for rcp_index in xrange(len(RCP_list)):
        rcp = RCP_list[rcp_index]
        
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
        K_C_Reg = K_C_List[z]
        
        TempNS = concatenate((mean(TempNSdata[time<=2000+11/12])*ones(AddTime+12), TempNSdata)) # use mean(TempNSdata[time<=2000]) = mean WOD13 before year 2000 as initial temperature profile for spine up
        TempCORListCenter = (TempList - T0)/rho
        NormCor = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)
        Dist = norm.pdf(TempCORListCenter)*norm.cdf(TempCORListCenter*skew)/max(NormCor)
        maxD = max(Dist)
        T_opt = TempList[int(list(Dist).index(maxD))]  # Find the optimal temperature for coral growth (used to simulate occurence of bleaching)
        
        # uncomment 3 lines bellow to create the result data file, here I use a different variable for each filenames 
        
        fileCoral = open(folder+"CORAL-"+rcp+"-"+Location_value+".dat", "wr")
        fileTrait = open(folder+"TRAIT-"+rcp+"-"+Location_value+".dat", "wr")
        fileSymb = open(folder+"SYMB-"+rcp+"-"+Location_value+".dat", "wr")
        
        
        if Location_value == "GBR":
            N_List = (scale)*N_values
            Initial = Init["GBR"]
        elif Location_value == "SEA":
            N_List = (scale)*N_values
            Initial = Init["SEA"]
        elif Location_value== "CAR":
            N_List = (scale)*N_values
            Initial = Init["CAR"]
        else:
            print "That location is not on my list"
        for i in xrange(len(N_List)):
            N = N_List[i]
            ode15s = ode(SystemForcing)
            ode15s.set_f_params(T0, rho, skew, N, NormCor, TempNS, K_C_Reg, NumTime)
            # Choosing an integrator choosing a solver that can deal stiff problems since the monthly temperature forcing is not smooth
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
                        print rcp, Location_value, i, ode15s.t, "reset"
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
                        print rcp, Location_value, i, ode15s.t, "no reset"
                elif Location_value == "SEA": # bleaching thresholds for the SEA
                    if TempNS[int(ode15s.t)] >= T_opt+1: # bleaching occurs symbionts are expeled and the system is reset to restart again
                        random_reduction = random.uniform(1-0.95, 1-0.25)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Location_value, i, ode15s.t, "reset"
                    else:
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Location_value, i, ode15s.t, "no reset"
                elif Location_value == "CAR": # bleaching thresholds for the CAR
                    if TempNS[int(ode15s.t)] >= T_opt+1 and TempNS[int(ode15s.t)] < T_opt+3: # bleaching occurs symbionts are expeled and the system is reset to restart again
                        random_reduction = random.uniform(1-0.95, 1-0.15)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Location_value, i, ode15s.t, "reset"
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
                        print rcp, Location_value, i, ode15s.t, "no reset"
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
        fileSymb.close()
        
        HostDic[rcp] = HOST
        TraitDic[rcp] = TRAIT
        SymbDic[rcp] = SYMB
    
    return HostDic, TraitDic, SymbDic
        


