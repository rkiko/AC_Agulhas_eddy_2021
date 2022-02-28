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

I compute the oxygen respiration rate at different densities and over time from the oxygen variation.


\
an11

I take the PARR (particle associated respiration rate) already computed and (i) I transform in in micromol/kg/day; 
(ii) I plot its time series as a function of the depth


\
an12

I calculate the anomaly of chlorophyll along the BGC Argo float trajectory. For each profile/point of the Argo float
I extract (ia) the mean chl value in the mixed layer (iia) the maximum chlorophyll value in the mixed layer (iiia) the mean chlorophyll
value extracted by satellite inside the eddy maximal velocity contour (iva) and the maximum chlorophyll
value extracted by satellite inside the eddy maximal velocity contour; these values represent the local/eddy chlorophyll
content. Then, we fix a fixed circle centered at the location of the eddy, with radius equal to n (n at least 3) times the eddy radius, and,
for this circle, we extract (ib) the mean chlorophyll value inside the circle, and (iib) the mean chlorophyll value inside the circle
with the exclusion of the region covered by the eddy. The two latter values are considered as the environment chlorophyll content.
Then, we make the difference between the 4*2=8 combinations of these quantities and we obtain 8 chlorophyll anomalies.

\
an12_plot

Script in order to plot the data computed and saved with an12


\
an13

As in an09, I compute the respiration rate obtained from the oxygen trend as a function of the depth.
Here I calculate also the PARR as a function of the depth and I plot them together

\
an14

It is the same script as an05, with two main differences: (i) I exclude the dates 
after the float left the eddy; (ii) I plot the eddy radius and eddy-float distance as well

\
an15

It is the same script as an08, with two main differences: (i) I exclude the dates 
after the float left the eddy; (ii) I plot the eddy radius and eddy-float distance as well

\
an16_plot

I plot the data calculated with an12, but I do only one figure, which should be the final one used
in the paper. I plot the eddy trajectory colored proportionally to the chl anomaly (calculated as 
the difference between the satellite chl inside and outside the eddy), I plot the 
BGC Argo trajectory with the same color, I plot the BGC Argo profiles outside the eddy with a cross, and I plot
4 eddy contour, colored proportionally to the chl inside the eddy.

\
an17

It is the same script as an14 (plot of the times series of chl, temp etc.), with two main differences: (i) I exclude 
the profiles in which the BGC Argo float was outside the eddy; (ii) I do not plot the eddy radius and eddy-float distance as well

\
an18

It is the same script as an15 (plot of the times series of flux, Mip etc.+ integrated time series of sPOC and POC),
with two main differences: (i) I exclude the profiles in which the BGC Argo float was outside the eddy; (ii) I do not
plot the eddy radius and eddy-float distance as well


\
an19

It is the same than an13, but I also calculate the bbp respiration rate. To do so, first I calculate the
bbp (as POC) consumption in time. Then, I convert it to micromol/kg/day.
I also plot the time series of the bbp for each different isopcynal, as in an09

an19b

It is the same as an19, but I use a shorter time window (until the 23 June instead of mid August)


\
an20

I use the data calculated in an12 to plot the time series of the chlorophyll within the eddy.
The purpose is to see whether there was a chl bloom prior to the deployment of the BGC Argo float.
This could explain the presence of a patch of intense MiP and MaP concentration at about 300m depths at the beginning 
of the Argo deployment.


\
an21

I plot the time series of the eddy radius and of the distance of the BGC 
Argo float from the eddy center


\
an22

It is the same as an19, but I smooth the bbp with the savgol filter before the interpolation

an22b

It is the same as an22, but I use a shorter time window (until the 23 June instead of mid August)

\
an23

I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth) and I save them

an23b

I use the data calculated in an23 to calculate the carbon budget

\
an24

I use the dataset of Mouw et al. (https://doi.org/10.1594/PANGAEA.855594) and I try to extract the POC flux in 
the AC agulhas eddy region in order to compare it with our POC flux

\
an25

As in an 23, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth). This time, instead than 
calculating them for one parameterisation only (i.e., one fixed layer, one fixed starting and ending dates),
I do a loop on different layers (basically every 100m between 200 and 600m) and I calculate
the theoretical carbon disappearing due to respiration (difference of delta flux and delta integrated POC)
and I compare it with the PARR and carbon consumption associated to oxygen consumption with a plot.

\
an26

I plot the diffPSD_slope_func


\
an27

It is exactly as an04, so that I plot the relationship between particle size and particle
sinking speed extracted from the profiles of particle abundance in time and depth.
The difference is that I use the data corrected for the particle concentration, 
so I use different thresholds. I also plot a function of the sinking speed vs size from an
a different work.


\
an28

I calculate, for different depth bins, the particle size distribution(PSD). Here, I consider 
all the profiles included between the 23 of April and the 18th of August (the final
date can be changed), so that the PSD is a mean. Finally, I extract the particle abundance for the
smaller size classes (for which the uvp does not work) from the fit.

\
an29

The principle of this script is the same as in an28, with the difference that the PSD slope is extracted 
for each profile separately. In this way, I can evaluate how the particle abundance for 
the smaller size classes changes with time (even if the uncertainty is larger)

\
an30

In this script I download the data from the transect between South Africa and South America,
I conver the raw data, I append the CTD data and I calculate the fluxes and respiration rates.

\
an31

In this script I take the flux data from the transect between South Africa and South America
and I compare it with the data from our cruise