import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import pickle
import cartopy
from datetime import date
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

### I load eddy radius
filename_eddyData='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies.csv' % home
data_eddy=pd.read_csv(filename_eddyData, sep=',', header=0)
radius_Vmax=data_eddy['Rmax']
del data_eddy

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
### I don't need to filter the eddy contour data as they are updated and contain only eddy contours of the period in which
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
sel_eddyCont_missing_days=sel_Satel2Float!=99999
lonEddy_4float=np.array(lonEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
latEddy_4float=np.array(latEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
lonVmax_4float=lonVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]
latVmax_4float=latVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]
radius_Vmax=radius_Vmax[sel_Satel2Float[sel_eddyCont_missing_days]]
chl_inside_mean_4float=chl_inside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]

# #######################################################################
# # I exclude the BGC Argo days for which I do not have an eddy contour (and thus satellite chlorophyll value)
# #######################################################################
#
# sel_eddyCont_missing_days=sel_Satel2Float!=99999
# chla_float_mean=chla_float_mean[sel_eddyCont_missing_days]
# chla_float_max=chla_float_max[sel_eddyCont_missing_days]
# lon_float=lon_float[sel_eddyCont_missing_days]
# lat_float=lat_float[sel_eddyCont_missing_days]
# Date_Num_float=Date_Num_float[sel_eddyCont_missing_days]
# Date_Vec_float=Date_Vec_float[sel_eddyCont_missing_days]
#
# chl_inside_mean_4float=chl_inside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]
# chl_inside_max_4float=chl_inside_max[sel_Satel2Float[sel_eddyCont_missing_days]]
#
# chl_outside_mean_4float=chl_outside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]
# chl_inside_and_outside_mean_4float=chl_inside_and_outside_mean[sel_Satel2Float[sel_eddyCont_missing_days]]
#
# lonEddy_4float=np.array(lonEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
# latEddy_4float=np.array(latEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
# lonVmax_4float=lonVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]
# latVmax_4float=latVmax[:,sel_Satel2Float[sel_eddyCont_missing_days]]

#######################################################################
# I calculate the distance between the BGC argo and the eddy center
#######################################################################
dist=np.arccos(np.sin(lat_float*np.pi/180)*np.sin(latEddy_4float*np.pi/180)+np.cos(lat_float*np.pi/180)*np.cos(latEddy_4float*np.pi/180)*np.cos((lonEddy_4float-lon_float)*np.pi/180))*180/np.pi
dist_km=dist*111

#######################################################################
# I exclude the profiles in which the BGC argo float was outside the eddy
#######################################################################
sel_insideEddy=dist_km<=radius_Vmax
lon_float_inside=lon_float[sel_insideEddy]
lat_float_inside=lat_float[sel_insideEddy]
lon_float_outside=lon_float[~sel_insideEddy]
lat_float_outside=lat_float[~sel_insideEddy]

#######################################################################################################################
############### I plot the BGC Argo float (F) trajectory with the associated anomaly. There are 2 values of local chlorophyll
############### and 2 of environmental chlorophyll, for a total of 4 anomalies.
#######################################################################################################################
def resize_colobar(event):
    plt.draw()

    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,0.04, posn.height])

anom_meanE_meanOut=chl_inside_mean-chl_outside_mean

set_xlim_lower, set_xlim_upper = 5,18.5#min(lonEddy.min(),lon_float.min())-window_margin1,max(lon_float.max(),lonEddy.max())+window_margin1
set_ylim_lower, set_ylim_upper = -37,-31#min(latEddy.min(),lat_float.min())-window_margin1,max(lat_float.max(),latEddy.max())+window_margin1
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5, label='Eddy trajectory')  # cmap='Blues_r')
ax.plot(0, 0, 'dk-',linewidth=1.5, label='BGC float trajectory')  # cmap='Blues_r')
plot1 = ax.scatter(lonEddy, latEddy, c=anom_meanE_meanOut, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5,zorder=0)  # cmap='Blues_r')
plot2 = ax.scatter(lon_float_inside, lat_float_inside, facecolor='none',marker='d',edgecolor='k',linewidth=1.5,s=65,zorder=5)  # cmap='Blues_r')
plot3 = ax.scatter(lon_float_outside, lat_float_outside, c='k',marker='x',zorder=4)  # cmap='Blues_r')
ax.plot(lon_float, lat_float, 'k', alpha=0.9,zorder=3)
plot2 = plt.scatter(0, 0, c=chl_inside_mean_4float[25],cmap='Greens', vmin=0, vmax=0.5, s=1)
a = plot2.to_rgba(chl_inside_mean_4float[25])
plt.plot(0, 0,c=a,linewidth=3, label='Eddy contour')
plt.scatter(0,0, facecolors=a, marker="^", s=100,edgecolors='k',linewidth=1.5, label='Eddy center')
for i in [10, 25, 38, 50]:
    plot2 = plt.scatter(lonVmax_4float[:, i], latVmax_4float[:, i],c=lonVmax_4float[:, i]*0+chl_inside_mean_4float[i], cmap='Greens',vmin=0,vmax=0.5,s=1)
    a=plot2.to_rgba(chl_inside_mean_4float[i])
    plt.plot(lonVmax_4float[:, i], latVmax_4float[:, i],c=a,linewidth=3,zorder=10)
    plt.scatter(lonEddy_4float[i], latEddy_4float[i], facecolors=a, marker="^", s=100,edgecolors='k',linewidth=1.5,zorder=20)
    plt.scatter(lon_float[i], lat_float[i], edgecolor=a, s=180, facecolors='none',linewidth=1.5,zorder=15)

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(plot1,cax=ax_cb)
ax_cb2 = divider.new_horizontal(size="5%", pad=1, axes_class=plt.Axes)
fig.add_axes(ax_cb2)
cbar2 = plt.colorbar(plot2,cax=ax_cb2)

cbar.ax.set_ylabel('Chlorophyll anomaly (mg/m$^3$)', fontsize=18)
cbar2.ax.set_ylabel('Chlorophyll inside eddy (mg/m$^3$)', fontsize=18)
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
gl.xlabel_style = {'fontsize': 15}
gl.ylabel_style = {'fontsize': 15}
ax.set_xlabel('Longitude', fontsize=15)
ax.set_ylabel('Latitude', fontsize=15)
gl.right_labels = False
gl.top_labels = False
ax.legend(fontsize=15)

fig.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an16/Eddy_Float_Trj_and_Chl_anom_an16.pdf' % (home))
plt.close()






