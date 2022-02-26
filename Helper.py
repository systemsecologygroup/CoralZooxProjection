##### Simple model runs ######
#from Model_RCP_ode_monthly_2 import* # Model run is without acclimation from year 2010
from Model_RCP_ode_monthly import *  # Model run normal with acclimation all years

import pdb

##### Whole ranges for the parameter N (Used to get Figure 2 in manuscript) 
#### Values are still scalled by 1e-11 in Model_RCP_ode_monthly or Model_RCP_ode_monthly_2
rawNum_GBR_orig = arange(0.01, 0.1, 0.00005)
rawNum_SEA_orig = arange(0.01, 0.1, 0.00005)
rawNum_CAR_orig = arange(0.01, 0.1, 0.00005) 

# Values used for simulation of dynamics in Results/ 
# change this part for other simulations ### 
# These indices are the ones producing the main results in the manuscript
N_indexes = {"GBR":908, "SEA":330, "CAR":275}

rawNum_GBR = array([rawNum_GBR_orig[N_indexes["GBR"]]]) 
rawNum_SEA = array([rawNum_SEA_orig[N_indexes["SEA"]]])
rawNum_CAR = array([rawNum_CAR_orig[N_indexes["CAR"]]])

N_values = {"GBR":rawNum_GBR, "SEA":rawNum_SEA, "CAR":rawNum_CAR}

RCP_list = ["RCP26", "RCP45", "RCP85"]

Locations = array(["GBR", "SEA", "CAR"])

# RUN_Simulations for each locations,
# the function RUN_SIM the result are saved in the folder Results/

for Locs in Locations:
    # Returns type dictionary with the name of RCP scenario in RCP_list as key for the corresponding scenario
    
    '''Model run is without acclimation from year 2010 --- deactivate acclimation from year 2010 (i.e. N = 0)'''
    #Coral_scenarios, Trait_scenarios, Symb_scenarios = RUN_SIM(RCP_list, Locs, N_values[Locs], folder = "Results-N-0-2010_2/") # change to appropriate folder
    
    '''Model run normal all years --- ATTENTION: change line 2 to line 3, change Model_RCP_ode_monthly_2 to Model_RCE_ode_monthly'''
    Coral_scenarios, Trait_scenarios, Symb_scenarios = RUN_SIM(RCP_list, Locs, N_values[Locs], folder = "Results_2/") # change to appropriate folder
    
    



