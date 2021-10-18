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

I compute the oxygen variation at a given density value in order to compute the respiration rate along this isopycnal.
I choose a start and end date in order to avoid considering the profiles in which the oxygen increases.

\
an07

In this script, I take the EcoPart data (i.e., the data downloaded from Ecopart throught the
paruvpy module functions, which contain information about the Flux, MiP and Map)
and the Coriolis data (i.e. downloaded through ftp, and which contain the BGC variables
such as temperature, chl-a etc.) and I create a unique data file which contain both the 
information about the MiP MaP etc. data, both the BGC data.

\
an08

This script is similar to an02, but it takes the new function in paruvpy which have been updated on 2021-10-15. It then 
plots the MiP, Map time series expressed as (i) number of particles (ii) mgC/m3, and the flux time series.
MiP, MaP, and Flux data are filtered with savgol function and then temporally and spatially interpolated

\
an09

I compute the oxygen variation at different density values in order to compute the respiration rate along these isopycnals.
I then transform each density value into a depth value by considering all the depths corresponding with this density
value for all the profiles, and by doing the mean value. 


\
an10

I compute the oxygen variation at different density and over time.
