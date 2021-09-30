# AC_Agulhas_eddy_2021
Project to conduct analysis of float data from WMO6903095, ensuse001b, SN P53337-20FR001, Vinci_Voyager_1, off SouthAfrica during SO283

----
Scripts description:

an01

First script to load Ecopart_mip_map_flux_data.tsv and plot the flux in the whole water column over time


an02

I do the scatterplot and the contourf of the flux interpolated. Before plotting, I filter the Flux values with a 
savgol_filter function. I do also the contourf of the MiP


an03

I plot the time series of NP 128 and 256 interpolated  and filtered with savgol function
