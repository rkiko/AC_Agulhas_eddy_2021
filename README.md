# AC_Agulhas_eddy_2021
Project to conduct analysis of float data from WMO6903095, ensuse001b, SN P53337-20FR001, Vinci_Voyager_1, off SouthAfrica during SO283

----
Scripts description:

an01

First script to load Ecopart_mip_map_flux_data.tsv and plot the flux in the whole water column over time

\
an02

I do the scatterplot and the contourf of the flux interpolated. Before plotting, I filter the Flux values with a 
savgol_filter function. I do also the contourf of the MiP

\
an03

I plot the time series of NP for several size classes interpolated and filtered with savgol function

\
an04

I take the time series of NP for several size classes and then, for each size class, I try to compute
the steepness of the abundance of particles in time, which should give us the settling velocity of that size class

\
an05

I analyse the file 6903095_Sprof.nc and I plot the time series of chla, oxygen, salinity, bbp, and temperature
(filtered with the savgol function and interpolated over time).

\
an06

I compute the oxygen variation at a given density in order to compute the respiration rate