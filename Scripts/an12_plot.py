import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import pickle
import cartopy
from pathlib import Path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#######################################################################################################################
############### Loading Data
#######################################################################################################################
a_file = open("%s/an12/data_an12.pkl" % storedir, "rb")
data_an12 = pickle.load(a_file)
chla_float_mean=data_an12['chla_float_mean']
chla_float_max=data_an12['chla_float_max']
lon_float=data_an12['lon_float']
lat_float=data_an12['lat_float']
Date_Num_float=data_an12['Date_Num_float']
Date_Vec_float=data_an12['Date_Vec_float']
chl_inside_mean=data_an12['chl_inside_mean']
chl_inside_max=data_an12['chl_inside_max']
chl_outside_mean=data_an12['chl_outside_mean']
chl_inside_and_outside_mean=data_an12['chl_inside_and_outside_mean']
Date_Num_Eddy=data_an12['Date_Num_Eddy']
DateTime_Eddy=data_an12['DateTime_Eddy']
lonEddy=data_an12['lonEddy']
latEddy=data_an12['latEddy']
a_file.close()

#######################################################################
# I select the satellite chlorophyll data in the day of the BGC Argo float profiles
#######################################################################

sel_Satel2Float=np.squeeze(np.zeros((1,Date_Num_float.size)))
i=0
for i in range(0,Date_Num_float.size):
    datetmp = np.squeeze(np.array([Date_Vec_float[i,0:3]]))
    datetmp = datetime.datetime(datetmp[0],datetmp[1],datetmp[2])
    idx=np.array(np.where(datetmp == DateTime_Eddy))
    if idx.size>1: print('error, i %d, ' % i, idx)
    if idx.size==0:
        print('Warning: missing eddy contour and center for %d-%d-%d' % (Date_Vec_float[i,0],Date_Vec_float[i,1],Date_Vec_float[i,2]))
        sel_Satel2Float[i] = 99999
        continue
    sel_Satel2Float[i] = idx[0]

sel_Satel2Float=sel_Satel2Float.astype(int)
#######################################################################
# I exclude the BGC Argo days for which I do not have an eddy contour (and thus satellite chlorophyll value)
#######################################################################

sel_eddyCont_missing_days=sel_Satel2Float!=99999
chla_float_mean=chla_float_mean[sel_eddyCont_missing_days]
chla_float_max=chla_float_max[sel_eddyCont_missing_days]
lon_float=lon_float[sel_eddyCont_missing_days]
lat_float=lat_float[sel_eddyCont_missing_days]
Date_Num_float=Date_Num_float[sel_eddyCont_missing_days]
Date_Vec_float=Date_Vec_float[sel_eddyCont_missing_days]

chl_inside_mean_4float=chl_inside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]
chl_inside_max_4float=chl_inside_max[sel_Satel2Float[sel_eddyCont_missing_days]]

chl_outside_mean_4float=chl_outside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]
chl_inside_and_outside_mean_4float=chl_inside_and_outside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]

#######################################################################################################################
############### I plot the BGC Argo float (F) trajectory with the associated anomaly. There are 2 valeus of local chlorophyll
############### and 2 of environmental chlorophyll, for a total of 4 anomalies.
#######################################################################################################################

anom_meanF_meanOut=chla_float_mean-chl_outside_mean_4float
anom_maxF_meanOut=chla_float_max-chl_outside_mean_4float
anom_meanF_meanAll=chla_float_mean-chl_inside_and_outside_mean_4float
anom_maxF_meanAll=chla_float_max-chl_inside_and_outside_mean_4float

# Parameters for the plot
width, height = 0.8, 0.7
window_margin1=5
set_xlim_lower, set_xlim_upper = lon_float.min()-window_margin1,lon_float.max()+window_margin1
set_ylim_lower, set_ylim_upper = lat_float.min()-window_margin1,lat_float.max()+window_margin1

################# Plot part
fig = plt.figure(1, figsize=(12,8))
#ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(set_xlim_lower, set_xlim_upper))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
#ax.coastlines('10m')
plot1=plt.scatter(lon_float,lat_float,c=anom_meanF_meanOut,cmap='RdBu_r')
cbar = plt.colorbar(plot1)
plt.title('Chl anomaly: mean chl BGC Argo vs chl outside eddy', fontsize=18)
cbar.ax.set_ylabel('Chlorophyll anomaly (mg/m$^3$)', fontsize=18)
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
gl.xlabel_style = {'fontsize': 18}
gl.ylabel_style = {'fontsize': 18}
ax.set_xlabel('Longitude', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
gl.xlabels_top = False
gl.ylabels_right = False
fig.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an12/Chl_anom01_mean_vs_outside_an12.pdf' % (home))
plt.close()

#######################################################################################################################
############### I plot the eddy (E) center trajectory with the associated anomaly. There are 2 valeus of local chlorophyll
############### and 2 of environmental chlorophyll, for a total of 4 anomalies.
#######################################################################################################################

anom_meanE_meanOut=chl_inside_mean-chl_outside_mean
anom_maxE_meanOut=chl_inside_max-chl_outside_mean
anom_meanE_meanAll=chl_inside_mean-chl_inside_and_outside_mean
anom_maxE_meanAll=chl_inside_max-chl_inside_and_outside_mean






