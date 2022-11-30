# AC_Agulhas_eddy_2021
Project to conduct analysis of float data from WMO6903095, ensuse001b, SN P53337-20FR001, Vinci_Voyager_1, off SouthAfrica during SO283

----
Scripts description:

\
Fig_Main_v01

Script in which I plot all the graphics of the main text. There is still no
Fig. 1 and Fig. 2a


\
Fig_Main_v02

I change Fig 2e, 3a, and 4 because I use no more the depths but the isopycnals 
to show the change in average POC, flux, and carbon budget


\
Fig_Main_v03

I change Fig 4 because I do not include anymore the extended carbon budget. Also,
I do not plot anymore the carbon budget vs density, but only vs depth. Instead, 
I plot the carbon budget of the second time window


\
Fig_Main_v04

I do not plot anymore the carbon budget of the second time window. I use the 
bbpPOC from Koestner et al. as well. I plot the carbon budget without 
including the bbpPOC as well to test for differences. I add the bbp PARR. I
add the uncertainty in the conversion from oxy cons rate to carb cons rate.
I fix a bug in the extended size spectrum calculation and plot


\
Fig_Main_v05

For the carbon budget (bulk POC removal rate) I don't use anymore the POC 
content at the beginning and the end of the time window considered, but I 
obtain this values from a linear fit


\
Fig_Main_v06

For the carbon budget (bulk POC removal rate) I don't use anymore the POC 
content at the beginning and the end of the time window considered, but I 
obtain this values from the average of the POC content over n days. I move to
the Supplementary the Fig. 4 calculated with Cetinic relationshipe (for bbp to
POC) and without bbpPOC. Conversely, I move back as main figure the Fig. 4 with
the extended size spectrum


\
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
plot the eddy radius and eddy-float distance as well\
Update of 03 05 2022: 2 main changes: 1) I use the bbp from the high 
resolution dataset and I remove the spikes larger than 100 mgC 2) I calculate
the POC budget after having interpolated the Mip Map and bbp. Before, for each 
profile, I was taking its mean value. However several profiles do not go deeper
than 400 m. In this way, taking the mean of that profile would imply to take the mean of
the 200-400 m, and not of the 200-600m. This makes the estimate biased, because
close to the surface the concentrations are higher. Interpolating allows us to take
into account this fact


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

As in an 25, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth).
I do a loop on different layers (basically every 100m between 200 and 600m) and I calculate
the theoretical carbon disappearing due to respiration (difference of delta flux and delta integrated POC)
and I compare it with the PARR and carbon consumption associated to oxygen consumption with a plot.
The difference with an25 script is that I use (i) the respiration and MiP and MaP budgets calcualted
considering also the smallest size classes; (ii) the flux calculated considering also the smallest size classes 
and a different eta and b coefficient

\
an32

In this script I take the flux data from the transect between South Africa and South America
and I compare it with the data from our cruise and with the independent trap sediments from 
Mouw et al. (see script an24). I calculate the flux for our transect and for the South
Africa South America transect in 3 ways (i) the old way, in which I do not consider smallest
size classes and old eta and b values; (ii) including smallest size classes, old eta and b values
(iii) including smallest size classes and new eta and b values

\
an33

It is the same script as an18 (plot of the times series of flux, Mip etc.+ integrated time series of sPOC and POC),
with the difference that I plot also these quantities calculated with (i) a 
different eta and b value; (ii) a different eta and b value and the smallest size classes
03-05-22: edited the script in order to include modifications done in an18
(calculation of the POC budget by interpolating and not profile by profile)

\
an34

As in an 31, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth).
I do a loop on different layers (basically every 100m between 200 and 600m) and I calculate
the theoretical carbon disappearing due to respiration (difference of delta flux and delta integrated POC)
and I compare it with the PARR and carbon consumption associated to oxygen consumption with a plot.
The difference with an31 script is that I consider two time windows: the first goes from the 
13 April to x, the second goes from x to the 23 of September, with x which is changed.

\
an34b

I go only untile the 9 of September rather than to the 23, to see if it improves
(it does not, so I do not use these plots in principle)

\
an35

I download the Euphotic layer depth data at 9km resolution at 8-day resolution
\
an35b

I download the Euphotic layer depth data at 9km resolution at day resolution

\
an36

I calculate the Euphotic layer depth along the float trajectory using the product
at 8-day resolution
\
an36b

I calculate the Euphotic layer depth along the float trajectory using the product
at daily resolution. I do not save the data cos there are too nans

\
an37

As in an 17, I plot of the times series of chl, temp,oxygen vs depth
with the difference that I also add the euphotic depth layer calculated in an36

\
an38

As in an33, I plot of the times series of flux, Mip etc.+ integrated time series of sPOC and POC),
with the difference that I also add the euphotic depth layer calculated in an36
03-05-22: edited the script in order to include modifications done in an18
(calculation of the POC budget by interpolating and not profile by profile)

\
an39

As in an 37, I plot of the times series of chl, temp,oxygen vs depth
with the difference that I also add the contour lines of the density

\
an40

I use the dataset from Bouman et al. (https://doi.org/10.5194/essd-10-251-2018)
to calculate the alpha^B in the region of our eddy

\
an41

I download the attenuation coefficient k490

\
an42

I calculate the attenuation coefficient k490 along the float trajectory using the product
at 8-day resolution downloaded in an42

\
an43

I take the surface irradiance downloaded from Eumetsat for every day between the 13
April and the 23 September and I reduce the size of the matrix by considering only 
the region in which the eddy was trapped, in an analogous way to what I did
for the attenuation coefficient k490 (an41) and the Euphotic layer depth (an35)

\
an44

I calculate the surface irradiance downloaded (per hour) along the float trajectory using the product
at daily resolution downloaded in an43 from Eumetsat


\
an45

I calculate the critical depth


\
an46

I plot the critical depth over the chlorphyll, temperature, doxy etc. 
profile vs depth and time (as for an17 and an37)


\
an47

I try to do some delayed-time quality control on lon lat and juld (ie time)


\
an48

I analyse data from Argo floats which were in the surrounding region of our 
float between the 13 April and the 5 April 2022, and which were validated
as delayed time quality. I plot temp vs pres vs salinity for each profile of our argo 
float, and then I superpose the profiles of the neighbor floats
The data were downloaded from https://dataselection.euro-argo.eu/


\
an49

I do the same thing then in an49, but only for the profiles which sampled below
1000m 


\
an50

I plot the time series of the flux at 200 and 600 m


\
an51

I calculate the sinking speed of the flux


\
an52

I analyse the CTD data token by the vessel in the days of the Argo float 
deployment


\
an53

As in an 34, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth).
I do a loop on different layers (basically every 100m between 200 and 600m) and I calculate
the theoretical carbon disappearing due to respiration (difference of delta flux and delta integrated POC)
and I compare it with the PARR and carbon consumption associated to oxygen consumption with a plot.
I consider two time windows: the first goes from the 
13 April to x, the second goes from x to the 23 of September, with x which is changed.
The difference with an34 is that here I use the bbp from the BGC Argo float (Coriolis)
and not from Ecopart, so that the bbp is in higher resolution


\
an54

I analyse the dissolved oxygen climatological data downloaded from MIMOC oxygen climatology
(https://cloud.geomar.de/s/qp4ddYBFxpH2yCa)


\
an55

I test the function to calculate the MiP and MaP flux separately


\
an56

I calculate the MiP and MaP flux at 200m separately


\
an57

I calculate the distance of the float from the eddy center


\
an58

I plot the temperature profiles vs the float distance from the eddy center


\
an59

I plot the temperature profiles vs the float distance from the eddy center as


\
an60

For different depth layers, I plot the temperature vs distance from the eddy 
center, with each dot colored proportionally to the date in which it was 
measured



\
an61

I plot the TS diagram for the profiles inside the cyclone, coloring them differently
according to the date in which they were carried

\
an62

As in an60, I plot the temperature vs distance from the eddy 
center, with each dot colored proportionally to the date in which it was 
measured. The difference is that, rather than doing it for different depth 
layers, I do it for different isopycnals

\
an63

As in an 53, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth).
The difference with an53 is that here I calculate the flux along the isopycnals 
rather than a fixed depth.


\
an64

matlab script in which I try to analyse the data Remi prepared on 
the 16 June 2022 in which he reanalised the eddy merging events, shape and 
trajectories.

\
an65

I remake figure 1 including two novel eddy trajectories (the two eddies which
merge with the main one). The first eddy is done by hand for the moment


\
an66

As in an18, I plot the UVP quantities (i.e. flux, Mip POC etc) vs depth 
and time. The difference is that I use the new data for the bbp plots.


\
an67

As in an63, I calculate all the metrics necessary to estimate the carbon budget (
the Flux vs depth and time, the integrated POC, the PARR vs depth).
The difference with an63 is that here I exclude the second time window cos 
it is useless, and I change the panel label


\
an68

I calculate the ML depth for each profile and the isopycnal that better defines
it, and I save the data


\
an69

I download the chlorophyll inside and outside the two eddies (one is the 
cyclone targetd dy the BGC Argo float, and the other is the eddy which 
merged with this cyclone between the 1—11 August 2021) 


\
an70

I plot Fig. 1: Eddy 1 and Eddy 2 trajectories with their chl anomaly and 
contours, and float trajectory. I also plot time series of satellite chl 
inside and outside the Eddy 1 and Eddy 2, and the time series of eddy radius
and float distance from centroid


\
an71

I calculate sink speed of spherical clay particles using Dioguardi function


