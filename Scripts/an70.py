import os,sys
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import pickle
import cartopy
import netCDF4 as nc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum

#######################################################################################################################
############### Loading Data
#######################################################################################################################

# I load BGC Argo float trajectory
filename='6903095_Sprof_all.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))
lon_float=np.array(ds.variables['LONGITUDE'])
lat_float=np.array(ds.variables['LATITUDE'])
Date_Num_float=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)

# I load eddy trajectories and chlorophyll inside and outside
a_file = open("%s/an69/data_an69.pkl" % storedir, "rb")
data_an69 = pickle.load(a_file)
chl_inside_mean1=data_an69['chl_inside_mean1']
chl_outside_mean1=data_an69['chl_outside_mean1']
chl_inside_and_outside_mean1=data_an69['chl_inside_and_outside_mean1']
Date_Num_Eddy1=data_an69['Date_Num_Eddy1']
lonEddy1=data_an69['lonEddy1']
latEddy1=data_an69['latEddy1']
chl_inside_mean2=data_an69['chl_inside_mean2']
chl_outside_mean2=data_an69['chl_outside_mean2']
chl_inside_and_outside_mean2=data_an69['chl_inside_and_outside_mean2']
Date_Num_Eddy2=data_an69['Date_Num_Eddy2']
lonEddy2=data_an69['lonEddy2']
latEddy2=data_an69['latEddy2']
a_file.close()

### I load eddy contour data
filename_xVMax_eddy1='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy1_lon_an64m.csv' % home
filename_yVMax_eddy1='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy1_lat_an64m.csv' % home
data_xVMax1=pd.read_csv(filename_xVMax_eddy1, sep=',', header=0)
data_yVMax1=pd.read_csv(filename_yVMax_eddy1, sep=',', header=0)
lonVmax1=data_xVMax1.values[:,:]
latVmax1=data_yVMax1.values[:,:]
filename_xVMax_eddy2='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy2_lon_an64m.csv' % home
filename_yVMax_eddy2='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy2_lat_an64m.csv' % home
data_xVMax2=pd.read_csv(filename_xVMax_eddy2, sep=',', header=0)
data_yVMax2=pd.read_csv(filename_yVMax_eddy2, sep=',', header=0)
lonVmax2=data_xVMax2.values[:,:]
latVmax2=data_yVMax2.values[:,:]

#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)
sel_insideEddy = data_dist_radius['sel_insideEddy']
sel_insideEddy = (Date_Num_float<=day_end_timeseries)&(sel_insideEddy==1)

lon_float_outside = lon_float[(Date_Num_float<=day_end_timeseries)&(~sel_insideEddy)]
lat_float_outside = lat_float[(Date_Num_float<=day_end_timeseries)&(~sel_insideEddy)]
lon_float = lon_float[Date_Num_float<=day_end_timeseries]
lat_float = lat_float[Date_Num_float<=day_end_timeseries]

sel_eddy1=Date_Num_Eddy1<day_end_timeseries
chl_inside_mean1=chl_inside_mean1[sel_eddy1]
chl_outside_mean1=chl_outside_mean1[sel_eddy1]
Date_Num_Eddy1=Date_Num_Eddy1[sel_eddy1]
lonEddy1=lonEddy1[sel_eddy1]
latEddy1=latEddy1[sel_eddy1]

# sel_eddy2=Date_Num_Eddy2<day_end_timeseries #it's all true so I don't select nothing for eddy2


# #######################################################################
# # I save the chl values for the latex document
# #######################################################################
# sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
# from write_latex_data import write_latex_data
# from matlab_datevec import matlab_datevec
#
# Date_Vec_Eddy=np.zeros((Date_Num_Eddy.size,6))
# for i in range(0,Date_Num_Eddy.size):
#     Date_Vec_Eddy[i,:] = matlab_datevec(Date_Num_Eddy[i])
#
# Date_Vec_Eddy = Date_Vec_Eddy.astype(int)
# filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
# argument = 'chl_eddy_20201018'
# arg_value=np.mean(chl_inside_mean[0:3])
# write_latex_data(filename,argument,'%0.2f' % arg_value)
# argument = 'chl_percentage_difference_inside_outside'
# arg_value = np.mean(chl_inside_mean[175:339]-chl_outside_mean[175:339])/np.mean(chl_outside_mean[175:339])*100
# write_latex_data(filename,argument,'%0.d' % arg_value)

#######################################################################################################################
############### I plot the float + eddy trajectories with chl content and anomalies
#######################################################################################################################

def resize_colobar(event):
    plt.draw()
    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,0.04, posn.height])

anom_meanE_meanOut1=chl_inside_mean1-chl_outside_mean1
anom_meanE_meanOut2=chl_inside_mean2-chl_outside_mean2
# anom_meanE_meanOut_4float=anom_meanE_meanOut[sel_Satel2Float[sel_eddyCont_missing_days]]
chl_max_plot=0.6
set_xlim_lower, set_xlim_upper = 4,18.5#min(lonEddy.min(),lon_float.min())-window_margin1,max(lon_float.max(),lonEddy.max())+window_margin1
set_ylim_lower, set_ylim_upper = -37,-30.5#min(latEddy.min(),lat_float.min())-window_margin1,max(lat_float.max(),latEddy.max())+window_margin1
# set_xlim_lower, set_xlim_upper = 5,22 # For the larger plot
# set_ylim_lower, set_ylim_upper = -38,-29 # For the larger plot
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5, label='Eddy 1 trajectory')  # cmap='Blues_r')
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='deepskyblue',linewidth=0.5, label='Eddy 2 trajectory')  # cmap='Blues_r')
ax.plot(0, 0, 'k-',linewidth=2.5, label='BGC float trajectory')  # cmap='Blues_r')
plot1 = ax.scatter(lonEddy1, latEddy1, c=anom_meanE_meanOut1, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5,zorder=2)  # cmap='Blues_r')
plot4 = ax.scatter(lonEddy2, latEddy2, c=anom_meanE_meanOut2, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='deepskyblue',linewidth=0.5,zorder=2)  # cmap='Blues_r')
plot3 = ax.scatter(lon_float_outside, lat_float_outside, c='k',marker='x', s=100,zorder=12)  # cmap='Blues_r')
ax.plot(lon_float, lat_float, 'k', alpha=0.9,zorder=15,linewidth=2.5)
plot2 = plt.scatter(0, 0, c=chl_inside_mean1[25],cmap='Greens', vmin=0, vmax=chl_max_plot, s=1)
a = plot2.to_rgba(chl_inside_mean1[25])
plt.plot(0, 0,c=a,linewidth=3, label='Eddy contours')
i=0
for i in [0,90,237,275,297]:#306
    x=lonVmax1[:,i];y=latVmax1[:,i]
    tmp1 = plt.scatter(x, y,c=lonVmax1[:, i]*0+chl_inside_mean1[i], cmap='Greens',vmin=0,vmax=chl_max_plot,s=1)
    a=tmp1.to_rgba(chl_inside_mean1[i])
    tmp2 = plt.plot(x,y,c=a,linewidth=3,zorder=10)
    if i in [297]:  tmp5 = plt.plot(x,y,'k--',linewidth=1,zorder=30,label='Eddy merge')
    tmp3 = plt.scatter(lonEddy1[i], latEddy1[i], c=anom_meanE_meanOut1[i], cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='k',linewidth=0.5,zorder=20)
    idtmp=np.where(y==np.nanmin(y))
    if i in [237]: idtmp = np.where(x == np.nanmax(x))
    if i in [275]: idtmp = np.where(y == np.nanmax(y))
    if i in [297]: idtmp = np.where(x == np.nanmin(x))
    x1=x[idtmp][0];y1=y[idtmp][0]
    tmp4 = plt.plot(np.array([lonEddy1[i],x1]), np.array([latEddy1[i],y1]), c='k',alpha=0.2, linewidth=2, zorder=4)

for i in [3,32,116,189,268,325]:
    x=lonVmax2[:,i];y=latVmax2[:,i]
    tmp1 = plt.scatter(x, y,c=lonVmax2[:, i]*0+chl_inside_mean2[i], cmap='Greens',vmin=0,vmax=chl_max_plot,s=1)
    a=tmp1.to_rgba(chl_inside_mean2[i])
    tmp2 = plt.plot(x,y,c=a,linewidth=3,zorder=9)
    tmp3 = plt.scatter(lonEddy2[i], latEddy2[i], c=anom_meanE_meanOut2[i], cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='blue',linewidth=0.5,zorder=20)
    idtmp=np.where(y==np.nanmax(y))
    if i in [116,325]: idtmp=np.where(y==np.nanmin(y))
    if i in [189]: idtmp=np.where(x==np.nanmin(x))
    x1=x[idtmp][0];y1=y[idtmp][0]
    tmp4 = plt.plot(np.array([lonEddy2[i],x1]), np.array([latEddy2[i],y1]), c='k',alpha=0.2, linewidth=2, zorder=4)

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
ax.legend(fontsize=15)#,ncol=2)

fig.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an70/Eddy_Float_Trjs_and_Chl_anom_an70.pdf' % (home))
plt.close()

#######################################################################
# I plot eddy 1 radius vs time and relative float position
#######################################################################

# I load eddy 1 radius
filename_eddy1Data='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv' % home
data_eddy1=pd.read_csv(filename_eddy1Data, sep=',', header=0)
Date_Num_Eddy1=data_eddy1['Datenum']   # mean radius_Vmax  = 60.25 km
radius_Vmax1=data_eddy1['Rad_max']   # mean radius_Vmax  = 60.25 km
# I load float distance from eddy1 centroid
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)
Date_Num_float = data_dist_radius['Datenum']
Distance_centroid = data_dist_radius['Distance_Centroid']
# I define the end of the time series
day_end_timeseries=np.array([2021,9,23])
day_end_timeseries=matlab_datenum(day_end_timeseries)
# I define the eddy merging period
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = 0,150
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_float[0], day_end_timeseries))
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
plt.plot(Date_Num_Eddy1,radius_Vmax1,'b',label='Eddy radius')
plt.scatter(Date_Num_float,Distance_centroid,c='k',label='Float—eddy centroid distance',marker='*')
plt.scatter(Date_Num_float[~sel_insideEddy],Distance_centroid[~sel_insideEddy],label='Profiles excluded',marker='o', facecolors='none', edgecolors='r',s=60)
# I set xticks
nxticks = 10
xticks = np.linspace(Date_Num_float[0], day_end_timeseries, nxticks)
xticklabels = []
for i in xticks:
    tmp=matlab_datevec(i).astype(int)
    xticklabels.append(datetime.datetime(tmp[0],tmp[1],tmp[2],tmp[4],tmp[4],tmp[5]).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# plt.legend(fontsize=12)
plt.ylabel('Radius (km)', fontsize=15)
ax.text(-0.075, 1.05, 'a', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
ax.legend(fontsize=15,ncol=3)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an70/EddyRadiusAndDist_vs_time_an70.pdf' ,dpi=200)
plt.close()






