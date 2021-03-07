# General description

A trait-based model of coral-algae symbiosis to study coral acclimation under different warming scenarios. The model is applied to three coral reef ecosystems of the Tropics: the Great Barrier Reef, the South East Asia, and the Caribbean. Coral growth is constrained by thermal tolerance curves and corals respond to changing environmental conditions via plastic phenotypic changes.

<p align="center">
  <img src="mapn.png" width="500">
</p>

# Coral acclimation

Acclimation is captured by assuming that the temporal dynamics of a physiological trait, reflecting the energy that corals invest in the symbiotic relationship, is proportional to the gradient of coral fitness (assumed equal to net coral growth). The constant of proportionality reflects the speed of coral acclimation, i.e. the speed with which corals move towards an optimal trait value, one that maximises fitness, under changing environmental temperature. 

# Technicalities
Python 2.7 is required

User can simply use `Helper.py` to run one or more simulations for varying speed of acclimation N. The results of the simulation are saved in folders `Results/`, change to your own folder names. `Plot-Main-Monthy.py` can be used to display the main result (Figure 3 in manuscript) 

The python code `Model_RCP_ode_monthly.py` produces times series of coral biomass, coral energy investment trait and symbiont biomass. 

The simulation takes as input several parameters (which are defined within code) and are forced by time-series of monthly temperature scenarios from input files of the type `.dat`.

The time-series for monthly temperature are included in this repository (`Monthly-SST-scenario.zip`)

The model includes the process of bleaching (i.e. the expulsion of algae by the corals when the enviromenta temperature hits a certain threshold), which is parameterised according to literature data (listed in `S-density-bleaching.pdf`).


# Relevant references

P. Abrams, H. Matsuda, and Y. Harada. [Evolutionarily unstable fitness maxima and stable fitness minima of continuous traits](https://link.springer.com/article/10.1007/BF01237642). *Evolutionary Ecology*, **7**:465–487, 1993.

A. Merico, J. Bruggeman, and K. Wirtz. [A trait-based approach for downscaling complexity in plankton ecosystem models](https://www.sciencedirect.com/science/article/pii/S0304380009003275). *Ecological Modelling*, **220**:3001–3010, 2009.

N. A. Raharinirina, G. Brandt, and A. Merico. [A trait-based model for describing the adaptive dynamics of coral-algae symbiosis](https://www.frontiersin.org/articles/10.3389/fevo.2017.00031/full). *Frontiers in Ecology and Evolution*, **5**:31, 2017.
