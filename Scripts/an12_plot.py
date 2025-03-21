import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import pickle
import cartopy
from datetime import date
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

### I load eddy contour data
filename_xVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Xcon_max.csv' % home
filename_yVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Ycon_max.csv' % home
data_xVMax=pd.read_csv(filename_xVMax, sep=',', header=None)
data_yVMax=pd.read_csv(filename_yVMax, sep=',', header=None)
lonVmax=data_xVMax.values[:,:]
latVmax=data_yVMax.values[:,:]


#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
day0_insideEddy=[2021,4,13]
dayf_insideEddy=[2021,9,27]
day0_insideEddy=date.toordinal(date(day0_insideEddy[0], day0_insideEddy[1], day0_insideEddy[2]))+366
dayf_insideEddy=date.toordinal(date(dayf_insideEddy[0], dayf_insideEddy[1], dayf_insideEddy[2]))+366

###BGC Argo data
sel_insideEddy = [False for i in range(lon_float.size)]
i=0
for i in range(0,lon_float.size):
    daytmp=date.toordinal(date(Date_Vec_float[i,0], Date_Vec_float[i,1], Date_Vec_float[i,2]))+366
    if (daytmp>=day0_insideEddy)&(daytmp<=dayf_insideEddy):
        sel_insideEddy[i] = True

chla_float_mean=chla_float_mean[sel_insideEddy]
chla_float_max=chla_float_max[sel_insideEddy]
lon_float=lon_float[sel_insideEddy]
lat_float=lat_float[sel_insideEddy]
Date_Num_float=Date_Num_float[sel_insideEddy]
Date_Vec_float=Date_Vec_float[sel_insideEddy,:]

###Satellite data: here the difference is that I consider also the days before the BGC argo entered the eddy
### I don't need to filter the eddy contour data as they are updated and contin only eddy contours of the period in which
### BGC Argo was inside the eddy
sel_insideEddy = [False for i in range(lonEddy.size)]
i=0
for i in range(0,lonEddy.size):
    if (Date_Num_Eddy[i]<=dayf_insideEddy):
        sel_insideEddy[i] = True

chl_inside_mean=chl_inside_mean[sel_insideEddy]
chl_inside_max=chl_inside_max[sel_insideEddy]
chl_outside_mean=chl_outside_mean[sel_insideEddy]
chl_inside_and_outside_mean=chl_inside_and_outside_mean[sel_insideEddy]
Date_Num_Eddy=Date_Num_Eddy[sel_insideEddy]
DateTime_Eddy=DateTime_Eddy[sel_insideEddy]
lonEddy=lonEddy[sel_insideEddy]
latEddy=latEddy[sel_insideEddy]

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

lonEddy_4float=np.array(lonEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
latEddy_4float=np.array(latEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
lonVmax_4float=lonVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]
latVmax_4float=latVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]

#######################################################################################################################
############### I plot the BGC Argo float (F) trajectory with the associated anomaly. There are 2 values of local chlorophyll
############### and 2 of environmental chlorophyll, for a total of 4 anomalies.
#######################################################################################################################

anom_meanF_meanOut=chla_float_mean-chl_outside_mean_4float
anom_maxF_meanOut=chla_float_max-chl_outside_mean_4float
anom_meanF_meanAll=chla_float_mean-chl_inside_and_outside_mean_4float
anom_maxF_meanAll=chla_float_max-chl_inside_and_outside_mean_4float

# Parameters for the plot
width, height = 0.8, 0.7
window_margin1=3
set_xlim_lower, set_xlim_upper = 5,20#min(lonEddy.min(),lon_float.min())-window_margin1,max(lon_float.max(),lonEddy.max())+window_margin1
set_ylim_lower, set_ylim_upper = -37,-31#min(latEddy.min(),lat_float.min())-window_margin1,max(lat_float.max(),latEddy.max())+window_margin1
titles_list=['mean chl BGC Argo vs chl outside eddy','max chl BGC Argo vs chl outside eddy','mean chl BGC Argo vs chl all region','max chl BGC Argo vs chl all region']
labels_list=['meanFloat_vs_outside','maxFloat_vs_outside','meanFloat_vs_allRegion','maxFloat_vs_allRegion']

################# Plot part
# I do the loop on the 4 plots
i_anom=0
for i_anom in range(0,titles_list.__len__()):
    if i_anom == 0: anom=anom_meanF_meanOut
    elif i_anom == 1:   anom=anom_maxF_meanOut
    elif i_anom == 2:   anom=anom_meanF_meanAll
    elif i_anom == 3:   anom=anom_maxF_meanAll

    fig = plt.figure(1, figsize=(12,8))
    #ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(set_xlim_lower, set_xlim_upper))
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
    ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
    #ax.coastlines('10m')
    plot1=plt.scatter(lon_float,lat_float,c=anom,cmap='RdBu_r',vmin=-0.2,vmax=0.2)#cmap='Blues_r')
    plt.plot(lon_float,lat_float,'k',alpha=0.2,label='BGC float trajectory')
    plt.plot(np.nan,np.nan,'k',label='Eddy max. vel. contour')
    plt.scatter(np.nan,np.nan, c='k', marker="^", s=20, label='Eddy center')
    ax.legend(fontsize=15)
    cbar = plt.colorbar(plot1)
    plt.title('Chl anomaly: %s' % titles_list[i_anom], fontsize=18)
    cbar.ax.set_ylabel('Chlorophyll anomaly (mg/m$^3$)', fontsize=18)
    gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
    gl.xlabel_style = {'fontsize': 18}
    gl.ylabel_style = {'fontsize': 18}
    ax.set_xlabel('Longitude', fontsize=18)
    plt.ylabel('Latitude', fontsize=18)
    gl.xlabels_top = False
    gl.ylabels_right = False
    for i in [10,25,38,50]:
        plt.plot(lonVmax_4float[:,i],latVmax_4float[:,i],'k')
        plt.scatter(lonEddy_4float[i],latEddy_4float[i],c='k',marker="^",s=20,label='Eddy center')
        plt.scatter(lon_float[i],lat_float[i],edgecolors='g',s=108, facecolors='none')

    fig.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an12/Chl_anom%02d_%s_an12.pdf' % (home,i_anom+1,labels_list[i_anom]))
    plt.close()

#######################################################################################################################
############### I plot the eddy (E) center trajectory with the associated anomaly. There are 2 values of local chlorophyll
############### and 2 of environmental chlorophyll, for a total of 4 anomalies.
#######################################################################################################################

anom_meanE_meanOut=chl_inside_mean-chl_outside_mean
anom_maxE_meanOut=chl_inside_max-chl_outside_mean
anom_meanE_meanAll=chl_inside_mean-chl_inside_and_outside_mean
anom_maxE_meanAll=chl_inside_max-chl_inside_and_outside_mean

# Parameters for the plot
titles_list=['mean chl inside Eddy (satellite) vs chl outside eddy','max chl inside Eddy (satellite) vs chl outside eddy','mean chl inside Eddy (satellite) vs chl all region','max chl inside Eddy (satellite) vs chl all region']
labels_list=['meanEddy_vs_outside','maxEddy_vs_outside','meanEddy_vs_allRegion','maxEddy_vs_allRegion']

################# Plot part
# I do the loop on the 4 plots
i_anom=0
for i_anom in range(0,titles_list.__len__()):
    if i_anom == 0: anom=anom_meanE_meanOut
    elif i_anom == 1:   anom=anom_maxE_meanOut
    elif i_anom == 2:   anom=anom_meanE_meanAll
    elif i_anom == 3:   anom=anom_maxE_meanAll

    fig = plt.figure(1, figsize=(12,8))
    #ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(set_xlim_lower, set_xlim_upper))
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
    ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
    #ax.coastlines('10m')
    plt.scatter(0,0,c=0.5,cmap='RdBu_r',vmin=-0.2,vmax=0.2,label='Eddy center positions')#cmap='Blues_r')
    plot1=plt.scatter(lonEddy,latEddy,c=anom,cmap='RdBu_r',vmin=-0.2,vmax=0.2)#cmap='Blues_r')
    ax.legend(fontsize=15)
    cbar = plt.colorbar(plot1)
    plt.title('Chl anomaly: %s' % titles_list[i_anom], fontsize=18)
    cbar.ax.set_ylabel('Chlorophyll anomaly (mg/m$^3$)', fontsize=18)
    gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
    gl.xlabel_style = {'fontsize': 18}
    gl.ylabel_style = {'fontsize': 18}
    ax.set_xlabel('Longitude', fontsize=18)
    plt.ylabel('Latitude', fontsize=18)
    gl.xlabels_top = False
    gl.ylabels_right = False
    fig.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an12/Chl_anom%02d_%s_an12.pdf' % (home,i_anom+5,labels_list[i_anom]))
    plt.close()





