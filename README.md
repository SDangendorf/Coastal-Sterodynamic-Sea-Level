# Coastal-Sterodynamic-Sea-Level
Contains scripts and data to run the EOF reocnstruction presented in Dangendorf et al. (2021), in review in Nature Climate Change
It requires the three data files: 
(1) the virtual stations (Virtual_Stations.mat) of coastal residual sea level (after removing barystatic, vertical land motion and barotropic atmospheric components) and barotropic atmospheric sea level
(2) the Steric Height (Steric-Height.mat)
(3) the corresponding time series at individual tide gauges (TG_Data.mat)

Also included are the final budget components (Budget_Terms.mat) for barystatic sea level (SIBSL), sterodynamic sea level from the EOF reconstruction (SISDSL), their sum (SIBud) and the observed sea level corrected for vertical land motion (SInovlm). All components are given as percentiles (2.5,50,97.5) of the 5000-member ensemble presented in the paper.  
