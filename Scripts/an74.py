
########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 01
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Fig.1
import numpy as np
import pandas as pd
import os,sys
import netCDF4 as nc
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename_coriolis='6903095_Sprof_all.nc'
######
import datetime
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable

#######################################################################
# Loading Data
#######################################################################

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
Date_Num_Eddy1=data_an69['Date_Num_Eddy1']
lonEddy1=data_an69['lonEddy1']
latEddy1=data_an69['latEddy1']
chl_inside_mean2=data_an69['chl_inside_mean2']
chl_outside_mean2=data_an69['chl_outside_mean2']
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

#######################################################################
# I plot the float + eddy trajectories with chl content and anomalies
#######################################################################

def resize_colobar(event):
    plt.draw()
    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,0.04, posn.height])

anom_meanE_meanOut1=chl_inside_mean1-chl_outside_mean1
anom_meanE_meanOut2=chl_inside_mean2-chl_outside_mean2
chl_max_plot=0.6
# set_xlim_lower, set_xlim_upper = 4,18.5
# set_ylim_lower, set_ylim_upper = -37,-30.5
set_xlim_lower, set_xlim_upper = 4,22 # For the larger plot
set_ylim_lower, set_ylim_upper = -38,-30.5 # For the larger plot
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5, label='Eddy trajectory')  # cmap='Blues_r')
ax.plot(0, 0, 'k-',linewidth=2.5, label='BGC float trajectory')  # cmap='Blues_r')
plot1 = ax.scatter(lonEddy1, latEddy1, c=anom_meanE_meanOut1, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5,zorder=2)  # cmap='Blues_r')
# plot4 = ax.scatter(lonEddy2, latEddy2, c=anom_meanE_meanOut2, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='goldenrod',linewidth=0.5,zorder=2)  # cmap='Blues_r')
plot3 = ax.scatter(lon_float_outside, lat_float_outside, c='red',marker='x', s=100,zorder=50)  # cmap='Blues_r')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9,edgecolor='gray')
ax.text(17.5,-31.5,'Cape\nColumbine',horizontalalignment='center',color='gray', bbox=props)
ax.plot([17.5,17.844],[-31.5,-32.8244],color='gray')
ax.text(20,-32.5,'South\nAfrica',color='white',fontsize=14.5,horizontalalignment='center')
# region
ax.plot(lon_float, lat_float, 'k', alpha=0.9,zorder=15,linewidth=2.5)
plot2 = plt.scatter(0, 0, c=chl_inside_mean1[25],cmap='Greens', vmin=0, vmax=chl_max_plot, s=1)
a = plot2.to_rgba(chl_inside_mean1[25])
i=0

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(plot1,cax=ax_cb)

cbar.ax.set_ylabel('Chlorophyll anomaly (mg/m$^3$)', fontsize=18)
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
gl.xlabel_style = {'fontsize': 15}
gl.ylabel_style = {'fontsize': 15}
ax.set_xlabel('Longitude', fontsize=15)
ax.set_ylabel('Latitude', fontsize=15)
gl.right_labels = False
gl.top_labels = False
ax.legend(fontsize=14.5,ncol=2,loc='lower right')

plt.savefig('../Plots/an74/Fig01_an74.pdf',dpi=200)
plt.close()
# endregion
# endregion
