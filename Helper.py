##### Simple model runs ######
from Model_RCP_ode_monthly import*

##### Whole ranges for the parameter N (Used to get Figure 2 in manuscript) 
#### Values are still rescalled by 1e-11
rawNum_GBR_orig = arange(0.01, 0.1, 0.00005)
rawNum_SEA_orig = arange(0.01, 0.1, 0.00005)
rawNum_CAR_orig = arange(0.01, 0.1, 0.00005) 

# Values used for simulation of dynamics in Results/ 
# change this part for other simulations ### 
# These indices are the ones producing the main results in the manuscript
N_indexes = {"GBR":920, "SEA":330, "CAR":275}

rawNum_GBR = array([rawNum_GBR_orig[N_indexes["GBR"]]])
rawNum_SEA = array([rawNum_SEA_orig[N_indexes["SEA"]]])
rawNum_CAR = array([rawNum_CAR_orig[N_indexes["CAR"]]])

N_values = {"GBR":rawNum_GBR, "SEA":rawNum_SEA, "CAR":rawNum_CAR}

RCP_list = ["RCP26", "RCP45", "RCP85"]

Locations = array(["GBR", "SEA", "CAR"])

# RUN_Simulations for each locations,
# the function RUN_SIM returns None but the result are saved in the folder Results/

for Locs in Locations:
    RUN_SIM(RCP_list, Locs, N_values[Locs], folder = "Results/")
    



