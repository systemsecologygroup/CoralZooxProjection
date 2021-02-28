# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean,var, ones, cumsum, random, sin, pi, load
from scipy.stats import norm
from numpy import amin, amax, meshgrid, arange, isnan, logical_not, interp, concatenate, arange, isnan
from scipy.integrate import ode
from matplotlib import pyplot as plt
import matplotlib
import numpy as np


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
def SystemForcing(t, y, T0, rho, skew, N, NormCor, TempNS, K_C_Reg):
    dSystem = zeros(len(y))
    coral = y[0]
    u = y[1]
    E = 1 - exp(-beta*u)
    symbiont = y[2]
    # Introducing temperature dependence
    #print t
    if t<=NumTime+AddTime:
        Temperature1 = TempNS[t] # this should rather be some function of t and not index but it works here because timestep = 1
    else:
        Temperature1 = TempNS[NumTime+AddTime] # there are index out of range sometimes, maybe because of the integrator, so I just make sure not to get those errors
    Tcenter = (Temperature1-T0)/rho
    Gx1Forcing = G_C*norm.pdf(Tcenter)*norm.cdf(Tcenter*skew)/max(NormCor) # divide by max(NormCor) to make sure that the maximum is G_C at T0
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

# Values used for simulation of dynamics in Results/ 
# I change this part for other simulations
rawNum_GBR = rawNum_GBR_orig
rawNum_SEA = rawNum_SEA_orig
rawNum_CAR = rawNum_CAR_orig

TempList = linspace(0, 75, 1500)
RCP_list = ["RCP26", "RCP45", "RCP85"]

Locations = array(["GBR", "SEA", "CAR"])

# Here I use a different variable for each filenames, but it doesn't matter otherwise. 
fileCoral = {}
fileTrait = {}
fileSymb = {}
file0 = {}
file1 = {}

for rcp_index in xrange(len(RCP_list)):
    rcp = RCP_list[rcp_index]
    if rcp_index == 0:
        z0 = 0
        z1 = len(Locations)
    elif rcp_index == 1:
        z0 = 0
        z1 = len(Locations)
    else :
        z0 = 0
        z1 = len(Locations)
    for z in xrange(z0, z1):
        sub1 = plt.subplot(3, 1, 1)
        sub2 = plt.subplot(3, 1, 2)
        sub3 = plt.subplot(3, 1, 3)
        #variable time
        file0[str(rcp_index)+str(z)] = open("Months-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        time = load(file0[str(rcp_index)+str(z)], allow_pickle = True)   
        file0[str(rcp_index)+str(z)].close()
        #Scenarios
        file1[str(rcp_index)+str(z)] = open("SST-"+Locations[z]+"-"+rcp+"-MPI"+".dat", "r")
        TempNSdata = load(file1[str(rcp_index)+str(z)], allow_pickle = True)
        file1[str(rcp_index)+str(z)].close()
        
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
        fileCoral[str(rcp_index)+str(z)] = open("Results/"+rcp+"/CORAL-"+rcp+"-"+Locations[z]+".dat", "wr")
        fileTrait[str(rcp_index)+str(z)] = open("Results/"+rcp+"/TRAIT-"+rcp+"-"+Locations[z]+".dat", "wr")
        fileSymb[str(rcp_index)+str(z)] = open("Results/"+rcp+"/SYMB-"+rcp+"-"+Locations[z]+".dat", "wr")
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
            # Choosing an integrator choosing a solver that can deal stiff problems since the monthly temperature forcing is not smooth
            ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
            ode15s.set_initial_value(Initial, 0)
            Dynamics = zeros((len(time)+AddTime+12, 3))  
            Dynamics[0, 0] = ode15s.y[0]
            Dynamics[0, 1] = ode15s.y[1]
            Dynamics[0, 2] = ode15s.y[2]
            k=1
            while ode15s.successful() and ode15s.t <= len(time)-2+AddTime+12:
                if z == 0:   # bleaching thresholds for the GBR
                    if TempNS[ode15s.t] >= T_opt+2 and TempNS[ode15s.t] < T_opt+4: # bleaching occurs symbionts are expeled and the system is reset to restart again
                        random_reduction = random.uniform(1-0.95, 1-0.10)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Locations[z], i, ode15s.t, "reset"
                    elif TempNS[ode15s.t] >= T_opt+4:
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
                        print rcp, Locations[z], i, ode15s.t, "no reset"
                elif z == 1: # bleaching thresholds for the SEA
                    if TempNS[ode15s.t] >= T_opt+1: # bleaching occurs symbionts are expeled and the system is reset to restart again
                        random_reduction = random.uniform(1-0.95, 1-0.25)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Locations[z], i, ode15s.t, "reset"
                    else:
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Locations[z], i, ode15s.t, "no reset"
                else: # bleaching thresholds for the CAR
                    if TempNS[ode15s.t] >= T_opt+1 and TempNS[ode15s.t] < T_opt+3: # bleaching occurs symbionts are expeled and the system is reset to restart again
                        random_reduction = random.uniform(1-0.95, 1-0.15)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                        print rcp, Locations[z], i, ode15s.t, "reset"
                    elif TempNS[ode15s.t] >= T_opt+3 and TempNS[ode15s.t] < T_opt+6:
                        random_reduction = random.uniform(1-0.95, 1-0.35)
                        ode15s.set_initial_value(array([ode15s.y[0], ode15s.y[1], (random_reduction)*ode15s.y[2]]), ode15s.t)
                        ode15s.integrate(ode15s.t+1)
                        Dynamics[k, 0] = ode15s.y[0]
                        Dynamics[k, 1] = ode15s.y[1]
                        Dynamics[k, 2] = ode15s.y[2]
                    elif TempNS[ode15s.t] >= T_opt+6:
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
                        print rcp, Locations[z], i, ode15s.t, "no reset"
                k+=1
            HOST.append(Dynamics[:, 0])
            TRAIT.append(Dynamics[:, 1])
            SYMB.append(Dynamics[:, 2])
        array(HOST).dump(fileCoral[str(rcp_index)+str(z)])
        array(TRAIT).dump(fileTrait[str(rcp_index)+str(z)])
        array(SYMB).dump(fileSymb[str(rcp_index)+str(z)])
        fileCoral[str(rcp_index)+str(z)].flush()
        fileTrait[str(rcp_index)+str(z)].flush()
        fileSymb[str(rcp_index)+str(z)].flush()
        fileCoral[str(rcp_index)+str(z)].close()
        fileTrait[str(rcp_index)+str(z)].close()
        fileSymb[str(rcp_index)+str(z)].close()



