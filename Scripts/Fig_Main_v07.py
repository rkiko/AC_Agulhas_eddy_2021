
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
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5, label='Eddy 1 trajectory')  # cmap='Blues_r')
ax.scatter(0, 0, c=0, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='goldenrod',linewidth=0.5, label='Eddy 2 trajectory')  # cmap='Blues_r')
ax.plot(0, 0, 'k-',linewidth=2.5, label='BGC float trajectory')  # cmap='Blues_r')
plot1 = ax.scatter(lonEddy1, latEddy1, c=anom_meanE_meanOut1, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='gray',linewidth=0.5,zorder=2)  # cmap='Blues_r')
plot4 = ax.scatter(lonEddy2, latEddy2, c=anom_meanE_meanOut2, cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='goldenrod',linewidth=0.5,zorder=2)  # cmap='Blues_r')
plot3 = ax.scatter(lon_float_outside, lat_float_outside, c='red',marker='x', s=100,zorder=50)  # cmap='Blues_r')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9,edgecolor='gray')
ax.text(17.5,-31.5,'Cape\nColumbine',horizontalalignment='center',color='gray', bbox=props)
ax.plot([17.5,17.844],[-31.5,-32.8244],color='gray')
ax.text(20,-32.5,'South\nAfrica',color='white',fontsize=14.5,horizontalalignment='center')
# region
ax.plot(lon_float, lat_float, 'k', alpha=0.9,zorder=15,linewidth=2.5)
plot2 = plt.scatter(0, 0, c=chl_inside_mean1[25],cmap='Greens', vmin=0, vmax=chl_max_plot, s=1)
a = plot2.to_rgba(chl_inside_mean1[25])
plt.plot(0, 0,c=a,linewidth=3, label='Eddy contours')
i=0
for i in [0,90,237,275,297]:
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
    tmp3 = plt.scatter(lonEddy2[i], latEddy2[i], c=anom_meanE_meanOut2[i], cmap='RdBu_r', vmin=-0.2, vmax=0.2,s=70,edgecolor='darkgoldenrod',linewidth=0.5,zorder=20)
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
ax.legend(fontsize=14.5,ncol=2,loc='lower right')

plt.savefig('../Plots/Fig_Main_v07/Fig01_v07.pdf',dpi=200)
plt.close()
# endregion
# endregion

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 02
########################################################################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
######### Fig. 02a
########################################################################################################################
# region Fig.2a
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

# I load eddy 1 radius
filename_eddy1Data='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv' % home
data_eddy1=pd.read_csv(filename_eddy1Data, sep=',', header=0)
Date_Num_Eddy1=data_eddy1['Datenum']   # mean radius_Vmax  = 60.25 km
radius_Vmax1=data_eddy1['Rad_max']   # mean radius_Vmax  = 60.25 km
# I load float distance from eddy1 centroid
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)
sel_insideEddy = data_dist_radius['sel_insideEddy']
Date_Num_float = data_dist_radius['Datenum']
Distance_centroid = data_dist_radius['Distance_Centroid']
# I define the end of the time series
day_end_timeseries=np.array([2021,9,23])
day_end_timeseries=matlab_datenum(day_end_timeseries)
sel_insideEddy = (Date_Num_float<=day_end_timeseries)&(sel_insideEddy==1)
# I define the eddy merging period
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)

# I plot
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
    xticklabels.append(datetime.datetime(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# plt.legend(fontsize=12)
plt.ylabel('Radius (km)', fontsize=15)
ax.text(-0.075, 1.05, 'a', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
ax.legend(fontsize=15,ncol=3)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig02a_v07.pdf',dpi=200)
plt.close()

#######################################################################
# I save some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
argument = 'dist_float_centroid_mean'
arg_value=np.mean(Distance_centroid[sel_insideEddy])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'mean_radius_0413_0731'
i1=np.where(Date_Num_Eddy1==matlab_datenum(2021,4,13))[0][0];i2=np.where(Date_Num_Eddy1==matlab_datenum(2021,7,31))[0][0]
radius_mean=np.mean(radius_Vmax1[i1:i2+1]);radius_std=np.std(radius_Vmax1[i1:i2+1])
write_latex_data(filename,argument,'%0.1f' % radius_mean)
argument = 'mean_surf_eddycore_0413_0731_km2'
arg_value=np.pi*(radius_mean*1000)**2/10**6
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'std_surf_eddycore_0413_0731_km2'
arg_value=2*np.pi*radius_mean*radius_std#np.pi*(radius_std*1000)**2/10**6
write_latex_data(filename,argument,'%d' % arg_value)

# endregion




########################################################################################################################
######### Fig. 02b,c,d,f
########################################################################################################################
# region Fig. 02b,c,d,f
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
import seawater as sw
import gsw
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to finish the time series in the plot
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)

#######################################################################
# I load the Coriolis data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename_coriolis))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")

#Standard variables
temp=np.array(ds.variables['TEMP_ADJUSTED'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])

#BGC Variables
chla=np.array(ds.variables['CHLA_ADJUSTED'])
doxy=np.array(ds.variables['DOXY_ADJUSTED'])
bbp700=np.array(ds.variables['BBP700_ADJUSTED'])

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
if np.sum(chla==99999)==chla.size:
    print('Taking non adjusted chlorophyll-a')
    chla = np.array(ds.variables['CHLA'])
if np.sum(doxy==99999)==doxy.size:
    print('Taking non adjusted oxygen')
    doxy = np.array(ds.variables['DOXY'])
if np.sum(bbp700==99999)==bbp700.size:
    print('Taking non adjusted bbp700')
    bbp700 = np.array(ds.variables['BBP700'])

#######################################################################
#I tranform the pressure to depth
#######################################################################
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth=np.ones(temp.shape)*99999
depth[mask_depth]=depth_tmp

#I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
#salinity and temperature first
mask_dens=np.logical_and(pres!=99999,temp!=99999,psal!=99999) # I exclude the points with value = 99999
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lon_tmp=np.tile(lon,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_dens]
lon_tmp=lon_tmp[mask_dens]
pres_tmp=pres[mask_dens]
psal_tmp=psal[mask_dens]
temp_tmp=temp[mask_dens]
abs_psal_tmp = gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp)  # I compute absolute salinity
cons_tmp = gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)          # I compute conservative temperature
dens_tmp = gsw.density.sigma0(abs_psal_tmp, cons_tmp)
abs_psal=np.ones(temp.shape)*99999
abs_psal[mask_dens]=abs_psal_tmp
cons_temp=np.ones(temp.shape)*99999
cons_temp[mask_dens]=cons_tmp
dens=np.ones(temp.shape)*99999
dens[mask_dens]=dens_tmp+1000

#######################################################################
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from oceanpy import bbp700toPOC
from oceanpy import bbp700toPOC_Koestner
sPOC=bbp700.copy()*0+99999
sPOC_Koestner=bbp700.copy()*0+99999
i=0
for i in range(0,bbp700.shape[0]):
    bbp700tmp=bbp700[i,:]
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
    chl_tmp=chla[i,:]
    # I exclude nan values
    sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)&(chl_tmp!=99999)
    bbp700tmp=bbp700tmp[sel]
    depth_tmp=depth_tmp[sel]
    temp_tmp=temp_tmp[sel]
    chl_tmp=chl_tmp[sel]
    # I convert to small POC (sPOC) and I set to 0 values <0
    sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
    sPOC_tmp[sPOC_tmp<0]=0
    sPOC[i,sel]=sPOC_tmp
    sPOC_tmp = bbp700toPOC_Koestner(bbp700tmp, chl_tmp)
    sPOC_tmp[np.isnan(sPOC_tmp)]=0
    sPOC_Koestner[i,sel]=sPOC_tmp

#######################################################################
# I select the data only when the BGC Argo float was inside the eddy AND before day_end_timeseries (which fixes the x limit)
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

lon=lon[sel_insideEddy]
lat=lat[sel_insideEddy]
Date_Num=Date_Num[sel_insideEddy]
pres=pres[sel_insideEddy]
depth=depth[sel_insideEddy,:]
dens=dens[sel_insideEddy,:]
cons_temp=cons_temp[sel_insideEddy,:]
abs_psal=abs_psal[sel_insideEddy,:]
chla=chla[sel_insideEddy,:]
doxy=doxy[sel_insideEddy,:]
sPOC=sPOC[sel_insideEddy,:]
sPOC_Koestner=sPOC_Koestner[sel_insideEddy,:]

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from oceanpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,chla.shape[0]):
    depth_tmp=depth[i,:]
    temp_tmp=cons_temp[i,:]
    # I exclude nan values
    sel_non_nan=(depth_tmp!=99999)&(temp_tmp!=99999)
    temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
    mld_tmp,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
    mld=np.append(mld,mld_tmp)

#######################################################################
# I load the critical depth
#######################################################################
a_file = open("%s/an45/data_an45.pkl" % storedir, "rb")
data_an45 = pickle.load(a_file)
critical_depth=data_an45['critical_depth']
critical_depth_datenum=data_an45['critical_depth_datenum']
critical_depth_datenum = critical_depth_datenum[~np.isnan(critical_depth)]
critical_depth = critical_depth[~np.isnan(critical_depth)]

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)

parameter_ylabel_list=['Temperature ($^{\circ}$C)','Chlorophyll $a$ (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)','$b_{bp}$POC (mgC $m^{-3}$)','Absolute salinity (g/kg)']
parameter_panellabel_list=['b','d','c','f','f','',]
parameter_savelabel_list=['b','d','c','f_Cetinic','f_Koestner','']
ipar=0
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0:   parameter=cons_temp.copy()
    elif ipar == 1: parameter=chla.copy()
    elif ipar == 2: parameter=doxy.copy()
    elif ipar == 3: parameter=sPOC.copy()
    elif ipar == 4: parameter=sPOC_Koestner.copy()
    elif ipar == 5: parameter=abs_psal.copy()

    #I filter the profiles
    parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        z = parameter[i, :]
        sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
        if ipar == 3: sel = (sel) & (z <= 100)
        if sum(sel) > 0:
            z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
            z = savgol_filter(z, 5, 1)
            parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
            depth_parameter = np.concatenate((depth_parameter, y1))
            dens_parameter = np.concatenate((dens_parameter, y2))

    parameter_filtered[parameter_filtered<0]=0
    dens_parameter[dens_parameter<1026]=1026
    dens_parameter[dens_parameter>1027.5]=1027.5
    # I define the x and y arrays for the contourf plot
    x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
    y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
    # I interpolate
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
    parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
    dens_interp_depth = griddata((Date_Num_parameter,depth_parameter), dens_parameter, (x_parameter_g, y_parameter_g), method="nearest")

    ########################################################
    ####### I plot: versus depth
    ########################################################
    if ipar==3:
        parameter_interp_depth[parameter_interp_depth > 40] = 40
    if ipar==4:
        parameter_interp_depth[parameter_interp_depth > 120] = 120

    width, height = 0.8, 0.7
    if ipar==0:
        set_ylim_lower, set_ylim_upper = y1_parameter.min(),700
    else:
        set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
    if (ipar==0):
        plt.plot(Date_Num,mld,'w')
        plot3 = ax.contour(x_parameter, y1_parameter, dens_interp_depth, levels=[1026.82,1027.2397618090454],colors='white', linestyles='dashed', linewidths=1, zorder=10)  # ,cmap='Blues_r')
        # ax.clabel(plot3, inline=1, fontsize=10)
        fmt = {}
        strs = ['1026.82 kg/m$^3$', '1027.24 kg/m$^3$']
        for l, s in zip(plot3.levels, strs):
            fmt[l] = s

        ax.clabel(plot3, plot3.levels[::], inline=True, fmt=fmt, fontsize=10)
    elif ipar==1:
        plt.plot(critical_depth_datenum,critical_depth,'w')#;plt.plot(critical_depth_datenum,critical_depth,'w.')

    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=700,color='w',linestyles='dashed')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=700,color='w',linestyles='dashed')
    # draw colorbar
    cbar = plt.colorbar(plot2)
    cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
    plt.ylabel('Depth (m)', fontsize=18)
    #plt.title('%smm' % NP_sizeclass, fontsize=18)
    #I set xticks
    nxticks=10
    xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
    xticklabels=[]
    for i in xticks:
        date_time_obj = date_reference + datetime.timedelta(days=i)
        xticklabels.append(date_time_obj.strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90,fontsize=12)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    # I save
    if ipar<5:
        plt.savefig('../Plots/Fig_Main_v07/Fig02%s_v07.pdf' % parameter_savelabel_list[ipar],dpi=200)
    else:
        plt.savefig('../Plots/Fig_Main_v07/Supplementary/SupplFig_salinity_time_depth_v07.pdf',dpi=200)
    plt.close()

# endregion





########################################################################################################################
######### Fig. 02e
########################################################################################################################
# region Fig. 02e

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
###########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to start and end the time series in the plot
day_start_timeseries=np.array([2021,4,13])
day_start_timeseries=matlab_datenum(day_start_timeseries)
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
ndays=int(day_end_timeseries-day_start_timeseries+1)

#######################################################################
# I load the mixed layer depth
#######################################################################
a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an68/data_MLD_an68.pkl" % (home), "rb")
data_an68 = pickle.load(a_file)
mld = data_an68["mld"]#"Date_Num": Date_Num,"lon": lon,"lat": lat}
mld_datenum = data_an68['Date_Num']
a_file.close()
#I convert datenum to calendar
i=0
for i in range(0,mld_datenum.size):
    date_time_obj = matlab_datevec(mld_datenum[i]).astype(int)
    date_time_obj = datetime.datetime(date_time_obj[0],date_time_obj[1],date_time_obj[2],date_time_obj[3],date_time_obj[4],date_time_obj[5])
    mld_datenum[i] = calendar.timegm(date_time_obj.timetuple())

#######################################################################
# I load and process MiP and MaP data
#######################################################################

#I load the file with the flux and POC
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data
Date_Time=np.array(data['Date_Time'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])
dens=np.array(data['Potential density [kg/m3]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.unique(Date_Num)

#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[0:sel_insideEddy.__len__()]
list_dates=list_dates[sel_insideEddy[0:list_dates.size]]
mld=mld[sel_insideEddy]
mld_datenum=mld_datenum[sel_insideEddy]
n_profiles=list_dates.size

#######################################################################
# I load and process bbp data
#######################################################################

# I load the bbp data and I select only those in the prescribed period
storedir = '%s/GIT/AC_Agulhas_eddy_2021/Data' % home
a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
data_an18 = pickle.load(a_file)
bbp_POC = data_an18['bbp_POC']
bbp_POC_Koestner = data_an18['bbp_POC_Koestner']
Date_Num_bbp = data_an18['Date_Num_bbp']
Date_Vec_bbp = data_an18['Date_Vec_bbp']
depth_bbp = data_an18['depth_bbp']
dens_bbp = data_an18['dens_bbp']
a_file.close()

sel_dates = sel_insideEddy[0:Date_Num_bbp.size]
Date_Num_bbp = Date_Num_bbp[sel_dates]
Date_Vec_bbp = Date_Vec_bbp[sel_dates,:]
depth_bbp = depth_bbp[sel_dates,:]
dens_bbp = dens_bbp[sel_dates,:]
bbp_POC = bbp_POC[sel_dates, :]
bbp_POC_Koestner = bbp_POC_Koestner[sel_dates, :]

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num_bbp_calendar = Date_Num_bbp.copy()
for i in range(0,Date_Num_bbp_calendar.size):
    date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                             Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
    Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.utcfromtimestamp(Date_Num[i])

########################################################################################################################
# Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
# interpolate it over a regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because
# some profiles only reach 400 m; (iii) for each density value of the density grid, I calculate the corresponding depth.
# Finally, (iv) I extract the mean MiP+MaP+bbp values between depth0 and depthf and between day0 and dayf (I obtain a
# time series)
########################################################################################################################

##############################################
# Step 1 and 2, filter and interpolation
MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])
bbp_Koestner_filtered=np.array([]);dens_bbp_Koestner_filtered=np.array([]);Date_Num_bbp_Koestner_filtered=np.array([])
MiP_MLD=np.array([]);MaP_MLD=np.array([]);bbp_MLD=np.array([]);bbp_Koestner_MLD=np.array([])
MiP_outMLD=np.array([]);MaP_outMLD=np.array([]);bbp_outMLD=np.array([]);bbp_Koestner_outMLD=np.array([])
isopycnal_600m = 1027.2397618090454 #calculated at step 4 of Fig. 3a

i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=MiP_POC[sel];x=Date_Num[sel];y=dens[sel];y2=depth[sel];sel2=~np.isnan(z);z=z[sel2];x=x[sel2];y=y[sel2];y2=y2[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        MiP_filtered = np.concatenate((MiP_filtered, z))
        Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x))
        dens_MiP_filtered = np.concatenate((dens_MiP_filtered, y))
        sel_MLD = y2<=mld[i]
        if sum(sel_MLD) > 0:
            MiP_MLD=np.append(MiP_MLD,np.mean(z[sel_MLD]))
        sel_outMLD = (y2>mld[i])&(y<isopycnal_600m)
        if sum(sel_outMLD) > 0:
            MiP_outMLD=np.append(MiP_outMLD,np.mean(z[sel_outMLD]))

    z=MaP_POC[sel];x=Date_Num[sel];y=dens[sel];y2=depth[sel];sel2=~np.isnan(z);z=z[sel2];x=x[sel2];y=y[sel2];y2=y2[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        MaP_filtered = np.concatenate((MaP_filtered, z))
        Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x))
        dens_MaP_filtered = np.concatenate((dens_MaP_filtered, y))
        sel_MLD = y2<=mld[i]
        if sum(sel_MLD) > 0:
            MaP_MLD=np.append(MaP_MLD,np.mean(z[sel_MLD]))
        sel_outMLD = (y2>mld[i])&(y<isopycnal_600m)
        if sum(sel_outMLD) > 0:
            MaP_outMLD=np.append(MaP_outMLD,np.mean(z[sel_outMLD]))

i=0
for i in range(0, bbp_POC.shape[0]):
    #Cetinic
    z=bbp_POC[i,:];y=dens_bbp[i,:];y2=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
    z[z>100] = 99999
    sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y=y[sel2];y2=y2[sel2]
    sel3=z==0
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        z[sel3]=0
        bbp_filtered = np.concatenate((bbp_filtered, z))
        Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, np.tile(x,sum(sel2)) ))
        dens_bbp_filtered = np.concatenate((dens_bbp_filtered, y))
        sel_MLD = y2<=mld[i]
        if sum(sel_MLD) > 0:
            bbp_MLD=np.append(bbp_MLD,np.mean(z[sel_MLD]))
        sel_outMLD = (y2>mld[i])&(y<isopycnal_600m)
        if sum(sel_outMLD) > 0:
            bbp_outMLD=np.append(bbp_outMLD,np.mean(z[sel_outMLD]))
    #Koestner
    z=bbp_POC_Koestner[i,:];y=dens_bbp[i,:];y2=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
    z[z>100] = 99999
    sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y=y[sel2];y2=y2[sel2]
    sel3=z==0
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        z[sel3]=0
        bbp_Koestner_filtered = np.concatenate((bbp_Koestner_filtered, z))
        Date_Num_bbp_Koestner_filtered = np.concatenate((Date_Num_bbp_Koestner_filtered, np.tile(x,sum(sel2)) ))
        dens_bbp_Koestner_filtered = np.concatenate((dens_bbp_Koestner_filtered, y))
        sel_MLD = y2<=mld[i]
        if sum(sel_MLD) > 0:
            bbp_Koestner_MLD=np.append(bbp_Koestner_MLD,np.mean(z[sel_MLD]))
        sel_outMLD = (y2>mld[i])&(y<isopycnal_600m)
        if sum(sel_outMLD) > 0:
            bbp_Koestner_outMLD=np.append(bbp_Koestner_outMLD,np.mean(z[sel_outMLD]))


# I define the x and y arrays for the MiP+MaP+bbp interpolation
x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), 100)
y_filtered = np.linspace(dens_bbp_filtered.min(), dens_MaP_filtered.max(), 200)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
MiP_interp = griddata((Date_Num_MiP_filtered, dens_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
MaP_interp = griddata((Date_Num_MaP_filtered, dens_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
bbp_interp = griddata((Date_Num_bbp_filtered, dens_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")
bbp_Koestner_interp = griddata((Date_Num_bbp_Koestner_filtered, dens_bbp_Koestner_filtered), bbp_Koestner_filtered,(x_filtered_g, y_filtered_g), method="nearest")

##############################################
# Step 3: for each density layer I calculate the corresponding mean depth
sel = (Date_Num <= list_dates.max())&(Date_Num > list_dates.min())
dens_tmp=dens[sel];depth_tmp=depth[sel]
sel2=(~np.isnan(dens_tmp))&(~np.isnan(depth_tmp))&(dens_tmp!=99999)&(depth_tmp!=99999)
dens_tmp=dens_tmp[sel2];depth_tmp=depth_tmp[sel2]

sel_bbp=(~np.isnan(dens_bbp))&(~np.isnan(depth_bbp))&(dens_bbp!=99999)&(depth_bbp!=99999)
dens_tmp_bbp=dens_bbp[sel_bbp];depth_tmp_bbp=depth_bbp[sel_bbp]

depth_bbp_filtered=np.array([]);depth_MiP_filtered=np.array([])
i=0
for i in range(0,y_filtered.size):
    if i==0:    dens0=y_filtered[i]
    else:       dens0=abs( y_filtered[i]+y_filtered[i-1] )*0.5
    if i==(y_filtered.size-1):  dens1=y_filtered[i]
    else:                       dens1=abs( y_filtered[i]+y_filtered[i+1] )*0.5
    #Depth calculation for MiP and MaP
    sel_dens = (dens_tmp >= dens0) & (dens_tmp < dens1)
    if sum(sel_dens) > 0:
        depth_MiP_filtered = np.append(depth_MiP_filtered, np.nanmean(depth_tmp[sel_dens]))
    else:
        depth_MiP_filtered = np.append(depth_MiP_filtered, np.array([np.nan]))

    # Depth calculation for bbp
    sel_dens = (dens_tmp_bbp >= dens0) & (dens_tmp_bbp < dens1)
    if sum(sel_dens) > 0:
        depth_bbp_filtered = np.append(depth_bbp_filtered, np.nanmean(depth_tmp_bbp[sel_dens]))
    else:
        depth_bbp_filtered = np.append(depth_bbp_filtered, np.array([np.nan]))



##############################################
# Step 4, I calculate the mean MiP+MaP+bbp (and std) between depth0 and depthf between day0 and dayf
tmp = np.abs(y_filtered - 1026.35)
idx1 = np.where( tmp == tmp.min() )[0][0]+1
tmp = np.abs(y_filtered - isopycnal_600m)
idx2 = np.where( tmp == tmp.min() )[0][0]
tmp = np.abs(y_filtered - 1026.82)
idx3 = np.where( tmp == tmp.min() )[0][0]+1

MiP_POC_0_102635=np.mean(MiP_interp[0:idx1,:],0)
MaP_POC_0_102635=np.mean(MaP_interp[0:idx1,:],0)
bbp_POC_0_102635=np.mean(bbp_interp[0:idx1,:],0)
bbp_POC_Koestner_0_102635=np.mean(bbp_Koestner_interp[0:idx1,:],0)

MiP_POC_0_102635_std = np.std(MiP_interp[0:idx1, :], 0)
MaP_POC_0_102635_std = np.std(MaP_interp[0:idx1, :], 0)
bbp_POC_0_102635_std = np.std(bbp_interp[0:idx1, :], 0)
bbp_POC_Koestner_0_102635_std = np.std(bbp_Koestner_interp[0:idx1, :], 0)

MiP_POC_102635_600=np.mean(MiP_interp[idx1:idx2,:],0)
MaP_POC_102635_600=np.mean(MaP_interp[idx1:idx2,:],0)
bbp_POC_102635_600=np.mean(bbp_interp[idx1:idx2,:],0)
bbp_POC_Koestner_102635_600=np.mean(bbp_Koestner_interp[idx1:idx2,:],0)

MiP_POC_102635_600_std = np.std(MiP_interp[idx1:idx2, :], 0)
MaP_POC_102635_600_std = np.std(MaP_interp[idx1:idx2, :], 0)
bbp_POC_102635_600_std = np.std(bbp_interp[idx1:idx2, :], 0)
bbp_POC_Koestner_102635_600_std = np.std(bbp_Koestner_interp[idx1:idx2, :], 0)

MiP_POC_102682_600=np.mean(MiP_interp[idx3:idx2,:],0)
MaP_POC_102682_600=np.mean(MaP_interp[idx3:idx2,:],0)
bbp_POC_102682_600=np.mean(bbp_interp[idx3:idx2,:],0)
bbp_POC_Koestner_102682_600=np.mean(bbp_Koestner_interp[idx3:idx2,:],0)

MiP_POC_102682_600_std = np.std(MiP_interp[idx3:idx2, :], 0)
MaP_POC_102682_600_std = np.std(MaP_interp[idx3:idx2, :], 0)
bbp_POC_102682_600_std = np.std(bbp_interp[idx3:idx2, :], 0)
bbp_POC_Koestner_102682_600_std = np.std(bbp_Koestner_interp[idx3:idx2, :], 0)

Integrated_POC_mgC_m3_0_102635 = MiP_POC_0_102635 + MaP_POC_0_102635 + bbp_POC_0_102635
Integrated_POC_mgC_m3_0_102635_std = np.sqrt( MiP_POC_0_102635_std**2 + MaP_POC_0_102635_std**2 + bbp_POC_0_102635_std**2 )
Integrated_POC_mgC_m3_102635_600 = MiP_POC_102635_600 + MaP_POC_102635_600 + bbp_POC_102635_600
Integrated_POC_mgC_m3_102635_600_std = np.sqrt( MiP_POC_102635_600_std**2 + MaP_POC_102635_600_std**2 + bbp_POC_102635_600_std**2 )
Integrated_POC_mgC_m3_102682_600 = MiP_POC_102682_600 + MaP_POC_102682_600 + bbp_POC_102682_600
Integrated_POC_mgC_m3_102682_600_std = np.sqrt( MiP_POC_102682_600_std**2 + MaP_POC_102682_600_std**2 + bbp_POC_102682_600_std**2 )
Integrated_POC_Koestner_mgC_m3_0_102635 = MiP_POC_0_102635 + MaP_POC_0_102635 + bbp_POC_Koestner_0_102635
Integrated_POC_Koestner_mgC_m3_0_102635_std = np.sqrt( MiP_POC_0_102635_std**2 + MaP_POC_0_102635_std**2 + bbp_POC_Koestner_0_102635_std**2 )
Integrated_POC_Koestner_mgC_m3_102635_600 = MiP_POC_102635_600 + MaP_POC_102635_600 + bbp_POC_Koestner_102635_600
Integrated_POC_Koestner_mgC_m3_102635_600_std = np.sqrt( MiP_POC_102635_600_std**2 + MaP_POC_102635_600_std**2 + bbp_POC_Koestner_102635_600_std**2 )
Integrated_POC_Koestner_mgC_m3_102682_600 = MiP_POC_102682_600 + MaP_POC_102682_600 + bbp_POC_Koestner_102682_600
Integrated_POC_Koestner_mgC_m3_102682_600_std = np.sqrt( MiP_POC_102682_600_std**2 + MaP_POC_102682_600_std**2 + bbp_POC_Koestner_102682_600_std**2 )
list_dates_Integrated_POC = x_filtered.copy()

# Integrated_POC_mgC_m3_0_MLD = MiP_MLD + MaP_MLD + bbp_MLD
# Integrated_POC_mgC_m3_MLD_600 = MiP_outMLD + MaP_outMLD + bbp_outMLD
# Integrated_POC_Koestner_mgC_m3_0_MLD = MiP_MLD + MaP_MLD + bbp_Koestner_MLD
# Integrated_POC_Koestner_mgC_m3_MLD_600 = MiP_outMLD + MaP_outMLD + bbp_Koestner_outMLD

# #######################################################################
# # I plot: Cetinic
# #######################################################################
# day_start_eddy_merging = datetime.datetime(2021,8,1)
# day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
# day_end_eddy_merging = datetime.datetime(2021,8,11)
# day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())
#
# width, height = 0.8, 0.5
# set_ylim_lower, set_ylim_upper = min(Integrated_POC_mgC_m3_0_102635.min(),Integrated_POC_mgC_m3_102682_600.min()*10),max(Integrated_POC_mgC_m3_0_102635.max()+Integrated_POC_mgC_m3_0_102635_std.max(),Integrated_POC_mgC_m3_102635_600.max()*10,Integrated_POC_mgC_m3_102682_600.max()*10)*1.02
#
# fig = plt.figure(1, figsize=(13,4))
# ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
# plt.plot(list_dates_Integrated_POC,Integrated_POC_mgC_m3_0_102635,'r',linewidth=3,label='0—1026.35 kg/m$^3$ [0—MLD]')
# plt.fill_between(list_dates_Integrated_POC, Integrated_POC_mgC_m3_0_102635 - Integrated_POC_mgC_m3_0_102635_std, Integrated_POC_mgC_m3_0_102635 + Integrated_POC_mgC_m3_0_102635_std,
#                   facecolor='r', color='r', alpha=0.2)#, label='Bulk POC\nremov. rate')
# # plt.plot(list_dates_Integrated_POC,Integrated_POC_mgC_m3_102635_600*10,'b',linewidth=3,label='1026.35—%0.2f kg/m$^3$ [MLD—600 m] ($\cdot$10)' % (isopycnal_600m))
# # plt.fill_between(list_dates_Integrated_POC, Integrated_POC_mgC_m3_102635_600*10 - Integrated_POC_mgC_m3_102635_600_std*0.5*10, Integrated_POC_mgC_m3_102635_600*10 + Integrated_POC_mgC_m3_102635_600_std*0.5*10,
# #                   facecolor='b', color='b', alpha=0.2)#, label='Bulk POC\nremov. rate')
# plt.plot(list_dates_Integrated_POC,Integrated_POC_mgC_m3_102682_600*10,'m',linewidth=3,label='1026.82—%0.2f kg/m$^3$ [Eddy Core] ($\cdot$10)' % (isopycnal_600m))
# plt.fill_between(list_dates_Integrated_POC, Integrated_POC_mgC_m3_102682_600*10 - Integrated_POC_mgC_m3_102682_600_std*0.5*10, Integrated_POC_mgC_m3_102682_600*10 + Integrated_POC_mgC_m3_102682_600_std*0.5*10,
#                   facecolor='m', color='m', alpha=0.2)#, label='Bulk POC\nremov. rate')
# plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# # I set xticks
# nxticks = 10
# xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
# xticklabels = []
# for i in xticks:
#     xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
# ax.set_xticks(xticks)
# ax.set_xticklabels(xticklabels)
# plt.xticks(rotation=90, fontsize=14)
# plt.legend(fontsize=12,ncol=2)
# plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
# ax.text(-0.075, 1.05, 'e', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
# plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# plt.savefig('../Plots/Fig_Main_v07/Fig02e_Cetinic_v07.pdf' ,dpi=200)
# plt.close()

#######################################################################
# I plot: Koestner
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(Integrated_POC_Koestner_mgC_m3_0_102635.min(),Integrated_POC_Koestner_mgC_m3_102682_600.min()*10),max(Integrated_POC_Koestner_mgC_m3_0_102635.max()+Integrated_POC_Koestner_mgC_m3_0_102635_std.max(),Integrated_POC_Koestner_mgC_m3_102682_600.max()*10)*1.02

fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates_Integrated_POC,Integrated_POC_Koestner_mgC_m3_0_102635,'r',linewidth=3,label='0—1026.35 kg/m$^3$ [0—MLD]')
plt.fill_between(list_dates_Integrated_POC, Integrated_POC_Koestner_mgC_m3_0_102635 - Integrated_POC_Koestner_mgC_m3_0_102635_std, Integrated_POC_Koestner_mgC_m3_0_102635 + Integrated_POC_mgC_m3_0_102635_std,
                  facecolor='r', color='r', alpha=0.2)#, label='Bulk POC\nremov. rate')
# plt.plot(list_dates_Integrated_POC,Integrated_POC_Koestner_mgC_m3_102635_600*10,'b',linewidth=3,label='1026.35—%0.2f kg/m$^3$ [MLD—600 m] ($\cdot$10)' % (isopycnal_600m))
# plt.fill_between(list_dates_Integrated_POC, Integrated_POC_Koestner_mgC_m3_102635_600*10 - Integrated_POC_Koestner_mgC_m3_102635_600_std*0.5*10, Integrated_POC_Koestner_mgC_m3_102635_600*10 + Integrated_POC_Koestner_mgC_m3_102635_600_std*0.5*10,
#                   facecolor='b', color='b', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.plot(list_dates_Integrated_POC,Integrated_POC_Koestner_mgC_m3_102682_600*10,'m',linewidth=3,label='1026.82—%0.2f kg/m$^3$ [Eddy Core] ($\cdot$10)' % (isopycnal_600m))
plt.fill_between(list_dates_Integrated_POC, Integrated_POC_Koestner_mgC_m3_102682_600*10 - Integrated_POC_Koestner_mgC_m3_102682_600_std*0.5*10, Integrated_POC_Koestner_mgC_m3_102682_600*10 + Integrated_POC_Koestner_mgC_m3_102682_600_std*0.5*10,
                  facecolor='m', color='m', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=12,ncol=2)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
ax.text(-0.075, 1.05, 'e', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig02e_Koestner_v07.pdf' ,dpi=200)
plt.close()

#######################################################################
# Supplementary Fig: time series of bbp, MiP, and MaP POC in the ML
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(MiP_POC_0_102635.min()*10,MaP_POC_0_102635.min()*10,bbp_POC_0_102635.min()),max(MiP_POC_0_102635.max()*10+MiP_POC_0_102635_std.max()*10,MaP_POC_0_102635.max()*10+MaP_POC_0_102635_std.max()*10,bbp_POC_Koestner_0_102635.max()+bbp_POC_Koestner_0_102635_std.max())*1.02

fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))

MiP_POC_0_102635 + MaP_POC_0_102635 + bbp_POC_Koestner_0_102635
plt.plot(list_dates_Integrated_POC,MiP_POC_0_102635*10,'r',linewidth=3,label='MiP POC ($\cdot$10)')
plt.fill_between(list_dates_Integrated_POC, MiP_POC_0_102635*10 - MiP_POC_0_102635_std*10, MiP_POC_0_102635*10 + MiP_POC_0_102635_std*10,
                  facecolor='r', color='r', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.plot(list_dates_Integrated_POC,MaP_POC_0_102635*10,'y',linewidth=3,label='MaP POC ($\cdot$10)')
plt.fill_between(list_dates_Integrated_POC, MaP_POC_0_102635*10 - MaP_POC_0_102635_std*10, MaP_POC_0_102635*10 + MaP_POC_0_102635_std*10,
                  facecolor='y', color='y', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.plot(list_dates_Integrated_POC,bbp_POC_Koestner_0_102635,'b',linewidth=3,label='$b_{bp}$POC')
plt.fill_between(list_dates_Integrated_POC, bbp_POC_Koestner_0_102635 - bbp_POC_Koestner_0_102635_std, bbp_POC_Koestner_0_102635 + bbp_POC_Koestner_0_102635_std,
                  facecolor='b', color='b', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=12,ncol=2)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
plt.title('Mixed Layer', fontsize=15, fontweight='bold')
ax.text(-0.075, 1.05, 'a', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/MiPMaPbbpPOC_Koestner_ML_v07.pdf' ,dpi=200)
plt.close()

#######################################################################
# Supplementary Fig: time series of bbp, MiP, and MaP POC in the Eddy Core
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(MiP_POC_102682_600.min()*1,MaP_POC_102682_600.min()*1,bbp_POC_102682_600.min()),20#max(MiP_POC_102682_600.max()*1+MiP_POC_102682_600_std.max()*1,MaP_POC_102682_600.max()*1+MaP_POC_102682_600_std.max()*1,bbp_POC_102682_600.max()+bbp_POC_102682_600_std.max())*1.02

fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))

MiP_POC_102682_600 + MaP_POC_102682_600 + bbp_POC_Koestner_102682_600
plt.plot(list_dates_Integrated_POC,MiP_POC_102682_600*1,'r',linewidth=3,label='MiP POC')
plt.fill_between(list_dates_Integrated_POC, MiP_POC_102682_600*1 - MiP_POC_102682_600_std*1, MiP_POC_102682_600*1 + MiP_POC_102682_600_std*1,
                  facecolor='r', color='r', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.plot(list_dates_Integrated_POC,MaP_POC_102682_600*1,'y',linewidth=3,label='MaP POC')
plt.fill_between(list_dates_Integrated_POC, MaP_POC_102682_600*1 - MaP_POC_102682_600_std*1, MaP_POC_102682_600*1 + MaP_POC_102682_600_std*1,
                  facecolor='y', color='y', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.plot(list_dates_Integrated_POC,bbp_POC_Koestner_102682_600,'b',linewidth=3,label='$b_{bp}$POC')
plt.fill_between(list_dates_Integrated_POC, bbp_POC_Koestner_102682_600 - bbp_POC_Koestner_102682_600_std, bbp_POC_Koestner_102682_600 + bbp_POC_Koestner_102682_600_std,
                  facecolor='b', color='b', alpha=0.2)#, label='Bulk POC\nremov. rate')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=12,ncol=2)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
plt.title('Eddy Core', fontsize=15, fontweight='bold')
ax.text(-0.075, 1.05, 'b', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/MiPMaPbbpPOC_Koestner_EC_v07.pdf' ,dpi=200)
plt.close()

#######################################################################
# I save values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=67;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'PercentageBbp_102682_102724_0413to0731'
arg_value=np.mean(bbp_POC_Koestner_102682_600[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_102682_600[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'PercentageBbp_0_102635_0413to0731'
arg_value=np.mean(bbp_POC_Koestner_0_102635[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_0_102635[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'PercentageMiP_0_102635_0413to0731'
arg_value=np.mean(MiP_POC_0_102635[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_0_102635[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'PercentageMiP_102682_600_0413to0731'
arg_value=np.mean(MiP_POC_102682_600[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_102682_600[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'PercentageMaP_0_102635_0413to0731'
arg_value=np.mean(MaP_POC_0_102635[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_0_102635[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'PercentageMaP_102682_600_0413to0731'
arg_value=np.mean(MaP_POC_102682_600[0:i])/np.mean(Integrated_POC_Koestner_mgC_m3_102682_600[0:i])*100
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'mean_Integrated_POC_Koestner_mgC_m3_102682_600_0413_0731'
mean_Integrated_POC_Koestner_mgC_m3_102682_600_0413_0731=np.mean(MiP_interp[idx3:idx2,0:i]) + np.mean(MaP_interp[idx3:idx2,0:i]) + np.mean(bbp_Koestner_interp[idx3:idx2,0:i])
write_latex_data(filename,argument,'%d' % mean_Integrated_POC_Koestner_mgC_m3_102682_600_0413_0731)
i=49;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_0_102635_0413to0703_value'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_0_102635[i])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'POC_0_102635_0413to0703_date'
write_latex_data(filename,argument,'3 July')
# Integrated_POC_mgC_m3_0_102635[49:58]
i=53;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_0_102635_10July'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_0_102635[i])
write_latex_data(filename,argument,'%0.1f' % arg_value)
i=67;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_0_102635_01August'
arg_value=np.mean([Integrated_POC_Koestner_mgC_m3_0_102635[i],Integrated_POC_Koestner_mgC_m3_0_102635[i-1]])
write_latex_data(filename,argument,'%0.1f' % arg_value)
i=73;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_0_102635_11August'
arg_value=np.mean([Integrated_POC_Koestner_mgC_m3_0_102635[i],Integrated_POC_Koestner_mgC_m3_0_102635[i-1]])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'POC_102682_102724_13April'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_102682_600[0])
write_latex_data(filename,argument,'%0.2f' % arg_value)

i=32;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
Integrated_POC_Koestner_mgC_m3_102682_600[26:34]
argument = 'POC_102682_102724_peak_date'
write_latex_data(filename,argument,'5 June')

Integrated_POC_Koestner_mgC_m3_102682_600[64:74]
i=69;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_102682_102724_date_start_increase'
write_latex_data(filename,argument,'5 August')
argument = 'POC_102682_102724_0812date'
i=73;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
arg_value=datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B')
write_latex_data(filename,argument,'%s' % arg_value)
argument = 'POC_102682_102724_0812value'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_102682_600[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)

i=60;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
Integrated_POC_Koestner_mgC_m3_102635_600[53:85]
argument = 'POC_102635_102724_0721date'
write_latex_data(filename,argument,'21 July')
argument = 'POC_102635_102724_0721value'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_102635_600[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=82;print(datetime.datetime.utcfromtimestamp(list_dates_Integrated_POC[i]).strftime('%d %B'))
argument = 'POC_102635_102724_0826_value'
arg_value=np.mean(Integrated_POC_Koestner_mgC_m3_102635_600[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)

#######Carbon subducted by the eddy
surfEddy_0413_0731_m2=18928197625.854225
subducting_speed_Eddy1_0413_0731_md = 1.37
Eddy_C_subrate = mean_Integrated_POC_Koestner_mgC_m3_102682_600_0413_0731*surfEddy_0413_0731_m2*subducting_speed_Eddy1_0413_0731_md/10**9
argument = 'Eddy_C_subrate_tons_perPers_perDay'
write_latex_data(filename,argument,'%0.1f' % Eddy_C_subrate)
CO2_emitted_EU2019 = 6.8/365
argument = 'nPeople_equivalentCO2emitted'
argvalue = Eddy_C_subrate/CO2_emitted_EU2019
write_latex_data(filename,argument,'%d' % argvalue)

plt.plot(x_filtered,bbp_POC_0_102635,'.b-')



# endregion


########################################################################################################################
######### Fig. 02g,h
########################################################################################################################
# region Fig. 02g,h
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
###########
import datetime
import seawater as sw
import gsw
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import calendar

#Here I define the time at which I want to finish the time series in the plot
# day_end_timeseries=datetime.datetime(2021,9,24)
# day_end_timeseries = calendar.timegm(day_end_timeseries.timetuple())
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)

#######################################################################
# I load the Coriolis data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename_coriolis))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num_bbp=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")

Date_Vec=np.zeros([Date_Num_bbp.size,6])
for i in range(0,Date_Num_bbp.size):
    date_time_obj = date_reference + datetime.timedelta(days=Date_Num_bbp[i]-matlab_datenum(1950,1,1))
    Date_Vec[i,0]=date_time_obj.year;Date_Vec[i,1]=date_time_obj.month;Date_Vec[i,2]=date_time_obj.day
    Date_Vec[i,3]=date_time_obj.hour;Date_Vec[i,4]=date_time_obj.minute;Date_Vec[i,5]=date_time_obj.second

Date_Vec=Date_Vec.astype(int)

#Standard variables
temp=np.array(ds.variables['TEMP_ADJUSTED'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])

#BGC Variables
chla=np.array(ds.variables['CHLA_ADJUSTED'])
doxy=np.array(ds.variables['DOXY_ADJUSTED'])
bbp700=np.array(ds.variables['BBP700_ADJUSTED'])

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
if np.sum(chla==99999)==chla.size:
    print('Taking non adjusted chlorophyll-a')
    chla = np.array(ds.variables['CHLA'])
if np.sum(doxy==99999)==doxy.size:
    print('Taking non adjusted oxygen')
    doxy = np.array(ds.variables['DOXY'])
if np.sum(bbp700==99999)==bbp700.size:
    print('Taking non adjusted bbp700')
    bbp700 = np.array(ds.variables['BBP700'])

#######################################################################
#I tranform the pressure to depth
#######################################################################
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth_bbp=np.ones(temp.shape)*99999
depth_bbp[mask_depth]=depth_tmp

#######################################################################
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from oceanpy import bbp700toPOC
bbp_POC=bbp700.copy()*0+99999
i=0
for i in range(0,bbp700.shape[0]):
    bbp700tmp=bbp700[i,:]
    depth_tmp=depth_bbp[i,:]
    temp_tmp=temp[i,:]
    # I exclude nan values
    sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)
    bbp700tmp=bbp700tmp[sel]
    depth_tmp=depth_tmp[sel]
    temp_tmp=temp_tmp[sel]
    # I convert to small POC (sPOC) and I set to 0 values <0
    sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
    sPOC_tmp[sPOC_tmp<0]=0
    bbp_POC[i,sel]=sPOC_tmp

#######################################################################
# I convert the bbp dates to float values (in seconds from 1970 1 1)
#######################################################################
Date_Num_bbp_calendar = Date_Num_bbp.copy()
for i in range(0, Date_Num_bbp_calendar.size):
    date_time_obj = datetime.datetime(Date_Vec[i, 0], Date_Vec[i, 1], Date_Vec[i, 2],
                             Date_Vec[i, 3], Date_Vec[i, 4], Date_Vec[i, 5])
    Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.utcfromtimestamp(Date_Num[i])

#######################################################################
# I load the MiP MaP data
#######################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data_ecopart=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data_ecopart.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data_ecopart
lon=np.array(data_ecopart['Longitude'][sel_filename])
lat=np.array(data_ecopart['Latitude'][sel_filename])
Date_Time=np.array(data_ecopart['Date_Time'][sel_filename])
pressure=np.array(data_ecopart['Pressure [dbar]'][sel_filename])
Flux=np.array(data_ecopart['Flux_mgC_m2'][sel_filename])
MiP_abund=np.array(data_ecopart['MiP_abun'][sel_filename])
MaP_abund=np.array(data_ecopart['MaP_abun'][sel_filename])
MiP_POC=np.array(data_ecopart['Mip_POC_cont_mgC_m3'][sel_filename])
MaP_POC=np.array(data_ecopart['Map_POC_cont_mgC_m3'][sel_filename])
depth=np.array(data_ecopart['Depth [m]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[sel_insideEddy[0:list_dates.size]]
Date_Num_bbp=Date_Num_bbp[sel_insideEddy]
Date_Num_bbp_calendar=Date_Num_bbp_calendar[sel_insideEddy]
depth_bbp=depth_bbp[sel_insideEddy]
temp=temp[sel_insideEddy]
bbp_POC=bbp_POC[sel_insideEddy,:]

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

ipar=0
parameter_shortname_list=['MiP_POC','MaP_POC','bbpPOC']
parameter_panellabel_list=['g','h','f']
parameter_ylabel_list=['MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)','$b_{bp}$POC (mgC $m^{-3}$)']
max_parameter_list=np.array([2.15,0.30,40])
MiP_POC_0_200=np.array([]);MiP_POC_200_600=np.array([])
MaP_POC_0_200=np.array([]);MaP_POC_200_600=np.array([])
bbp_POC_0_200=np.array([]);bbp_POC_200_600=np.array([])
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar == 0: parameter=MiP_POC.copy()
    elif ipar == 1: parameter=MaP_POC.copy()
    elif ipar == 2: parameter=bbp_POC.copy()

    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
    if ipar == 2:
        i=0
        for i in range(0, bbp_POC.shape[0]):
            z=parameter[i,:];y=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
            z[z>100] = 99999
            sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
            sel3 = z == 0
            if sum(sel2) > 0:
                z = savgol_filter(z, 5, 1)
                z[sel3] = 0
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, np.tile(x,sum(sel2)) ))
                depth_filtered = np.concatenate((depth_filtered, y2))
                # I define sel_200 and sel_200_600
                sel_0_200 = np.abs(y2) < 200
                sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
                bbp_POC_0_200=np.append(bbp_POC_0_200,np.mean(z[sel_0_200]));bbp_POC_200_600=np.append(bbp_POC_200_600,np.mean(z[sel_200_600]))
    else:
        # I filter the prophiles
        i=0
        for i in range(0,list_dates.size):
            sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
            z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
            if sum(sel2)>0:
                z=savgol_filter(z,5,1)
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
                depth_filtered = np.concatenate((depth_filtered, y2))
                # sel_200 and sel_200_600 are used only for the POC integrated in time
                sel_0_200 = np.abs(y2) < 200
                sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
                if ipar==0: MiP_POC_0_200=np.append(MiP_POC_0_200,np.mean(z[sel_0_200]));MiP_POC_200_600=np.append(MiP_POC_200_600,np.mean(z[sel_200_600]))
                if ipar==1: MaP_POC_0_200=np.append(MaP_POC_0_200,np.mean(z[sel_0_200]));MaP_POC_200_600=np.append(MaP_POC_200_600,np.mean(z[sel_200_600]))

    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    sel_0_200 = (np.abs(y_filtered) >= 0) & (np.abs(y_filtered) < 200)
    sel_200_600 = (np.abs(y_filtered) >= 200) & (np.abs(y_filtered) < 600)
    if ipar==0: MiP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MiP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==1: MaP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MaP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==2: bbp_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);bbp_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)

    if ipar == 2: continue

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
    # I draw colorbar
    cbar = plt.colorbar(plot2)
    cbar.ax.get_yticklabels()
    cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
    plt.ylabel('Depth (m)', fontsize=18)
    #I set xticks
    nxticks=10
    xticks=np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),nxticks)
    xticklabels=[]
    for i in xticks:
        xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90,fontsize=12)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/Fig_Main_v07/Fig02%s_v07.pdf' % (parameter_panellabel_list[ipar]),dpi=200)
    plt.close()


# endregion




########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 03
########################################################################################################################
########################################################################################################################
########################################################################################################################

#######################################################################################################################
######### Fig. 03a and Time Series of Flux in Eddy Core
########################################################################################################################
#region Fig. 03a and Time Series of Flux in Eddy Core
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
###########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
delta_depth_flux = 15  # around of the depth which I consider when extracting the flux

#Here I define the time at which I want to start and end the time series in the plot
day_start_timeseries=np.array([2021,4,13])
day_start_timeseries=matlab_datenum(day_start_timeseries)
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
ndays=int(day_end_timeseries-day_start_timeseries+1)

#######################################################################
# I load and process MiP and MaP data
#######################################################################

#I load the file with the flux and POC
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data
Date_Time=np.array(data['Date_Time'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])
dens=np.array(data['Potential density [kg/m3]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
Flux_MiP=np.array(data['Flux_MiP_mgC_m2'][sel_filename])
Flux_MaP=np.array(data['Flux_MaP_mgC_m2'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))

a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an68/data_MLD_an68.pkl" % (home), "rb")
data_an68 = pickle.load(a_file)
mld = data_an68['mld']
mld_datenum0 = data_an68['Date_Num']
a_file.close()
#I convert datenum to calendar
mld_datenum = mld_datenum0.copy()
i=0
for i in range(0,mld_datenum0.size):
    date_time_obj = matlab_datevec(mld_datenum0[i]).astype(int)
    date_time_obj = datetime.datetime(date_time_obj[0],date_time_obj[1],date_time_obj[2],date_time_obj[3],date_time_obj[4],date_time_obj[5])
    mld_datenum[i] = calendar.timegm(date_time_obj.timetuple())

#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[0:183]
list_dates=list_dates[sel_insideEddy[0:list_dates.size]]
mld=mld[sel_insideEddy]
mld_datenum=mld_datenum[sel_insideEddy]
mld_datenum0=mld_datenum0[sel_insideEddy]
n_profiles=list_dates.size

########################################################################################################################
# Here I extract the flux values at depth0 and depthf. To do so, (i) I filter it with a savgol function, then (ii) I
# interpolate it over a regular grid in time and density. This step is necessary to have the flux at 600 m, because some
# profiles only reach 400 m; (iii) for each density value of the density grid, I calculate the corresponding depth.
# Finally, (iv) I extract the flux values at depth0 and depthf
# I also extrac the flux, for each profile, at the corresponding MLD
########################################################################################################################

##############################################
# Step 1 and 2, filter and interpolation
Flux_filtered=np.array([]);Flux_MiP_filtered=np.array([]);Flux_MaP_filtered=np.array([]);dens_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([]);Flux_MLD=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=dens[sel];y2=depth[sel];z_MiP = Flux_MiP[sel];z_MaP = Flux_MaP[sel]
    sel2=~np.isnan(z);z=z[sel2];x=x[sel2];y=y[sel2];y2=y2[sel2];z_MiP=z_MiP[sel2];z_MaP=z_MaP[sel2]
    mld_tmp = mld[i]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        z_MiP = savgol_filter(z_MiP, 5, 1)
        z_MaP = savgol_filter(z_MaP, 5, 1)
        Flux_filtered = np.concatenate((Flux_filtered, z))
        Flux_MiP_filtered = np.concatenate((Flux_MiP_filtered, z_MiP))
        Flux_MaP_filtered = np.concatenate((Flux_MaP_filtered, z_MaP))
        Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x))
        dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y))
        sel_mld = (y2>=mld_tmp-delta_depth_flux)&(y2<mld_tmp+delta_depth_flux)
        if sum(sel_mld) > 0:
            Flux_MLD=np.append(Flux_MLD, np.mean(z[sel_mld]) )


# I define the x and y arrays for the Flux interpolation
x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
y_filtered = np.linspace(dens_Flux_filtered.min(), dens_Flux_filtered.max(), 200)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
Flux_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")
Flux_MiP_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
Flux_MaP_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")

##############################################
# Step 3: for each density layer I calculate the corresponding mean depth
depth_Flux_filtered=np.array([])
i=0
for i in range(0,y_filtered.size):
    if i==0:    dens0=y_filtered[i]
    else:       dens0=abs( y_filtered[i]+y_filtered[i-1] )*0.5
    if i==(y_filtered.size-1):  dens1=y_filtered[i]
    else:                       dens1=abs( y_filtered[i]+y_filtered[i+1] )*0.5
    #Depth calculation
    depth_Flux_filtered_tmp = np.array([])
    j=0
    for j in range(0,list_dates.size):
        sel=Date_Num==list_dates[j]
        dens_tmp=dens[sel];depth_tmp=depth[sel];sel2=(~np.isnan(dens_tmp))&(~np.isnan(depth_tmp));dens_tmp=dens_tmp[sel2];depth_tmp=depth_tmp[sel2]
        if sum(sel2) > 0:
            sel_dens = (dens_tmp >= dens0) & (dens_tmp < dens1)
            if sum(sel_dens) > 0:
                depth_Flux_filtered_tmp = np.append(depth_Flux_filtered_tmp, np.mean(depth_tmp[sel_dens]) )
            else:
                depth_Flux_filtered_tmp = np.append( depth_Flux_filtered_tmp, np.array([np.nan]) )
        else:
            depth_Flux_filtered_tmp = np.append( depth_Flux_filtered_tmp, np.array([np.nan]) )

    if sum(~np.isnan(depth_Flux_filtered_tmp))==0:
        depth_Flux_filtered = np.append(depth_Flux_filtered, np.array([np.nan]) )
    else:
        depth_Flux_filtered = np.append(depth_Flux_filtered, np.nanmean(depth_Flux_filtered_tmp))

##############################################
# Step 4, flux extraction at depth0 and depthf
depth0=200
depthf=600
sel_layer = (np.abs(depth_Flux_filtered) >= depth0-delta_depth_flux) & (np.abs(depth_Flux_filtered) < depth0+delta_depth_flux)
Flux_depth0 = np.mean(Flux_interp[sel_layer,:],axis=0)
Flux_MiP_depth0 = np.mean(Flux_MiP_interp[sel_layer,:],axis=0)
Flux_MaP_depth0 = np.mean(Flux_MaP_interp[sel_layer,:],axis=0)
isopycnal_depth0 = np.mean(y_filtered[sel_layer])
sel_layer = (np.abs(depth_Flux_filtered) >= depthf - delta_depth_flux) & (np.abs(depth_Flux_filtered) < depthf + delta_depth_flux)
Flux_depthf = np.mean(Flux_interp[sel_layer,:],axis=0)
Flux_MiP_depthf = np.mean(Flux_MiP_interp[sel_layer,:],axis=0)
Flux_MaP_depthf = np.mean(Flux_MaP_interp[sel_layer,:],axis=0)
isopycnal_depthf = np.mean(y_filtered[sel_layer])
depth_102635=np.interp(1026.35,y_filtered,depth_Flux_filtered)
sel_layer = (np.abs(depth_Flux_filtered) >= depth_102635 - delta_depth_flux) & (np.abs(depth_Flux_filtered) < depth_102635 + delta_depth_flux)
Flux_102635 = np.mean(Flux_interp[sel_layer,:],axis=0)

sel_layer= (y_filtered>isopycnal_depth0)&(y_filtered<=isopycnal_depthf)
Flux_eddy_core = np.mean(Flux_interp[sel_layer,:],axis=0)
#######################################################################
# I plot
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

fs=9
width, height = 0.82, 0.8

fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x_filtered,Flux_depth0,'r',label='%0.2f kg/m$^3$ [200 m]' % (isopycnal_depth0))
plt.plot(x_filtered,Flux_depthf,'b',label='%0.2f kg/m$^3$ [600 m]' % (isopycnal_depthf))
plt.plot(mld_datenum,Flux_MLD,'m',label='MLD')
# plt.plot(x_filtered,Flux_102635,'m',label='1026.35 kg/m$^3$')
plt.xlim(x_filtered.min(),x_filtered.max())
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
ax.text(-0.115, 1.05, 'a', transform=ax.transAxes,fontsize=14, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.ylabel('Flux(mgC/$m^2$/d)',fontsize=7)
plt.legend(fontsize=6,ncol=4)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
ax.set_xticks(xticks)
ax.set_xticklabels([])
plt.xticks(rotation=90, fontsize=7)
plt.savefig('../Plots/Fig_Main_v07/Fig03a_v07.pdf'  ,dpi=200)
plt.close()

# I transform the x_filtered in matlab datenum format
x_filtered_datenum = x_filtered.copy()
i = 0
for i in range(0, x_filtered.size):
    d = datetime.datetime.utcfromtimestamp(x_filtered[i])
    x_filtered_datenum[i] = matlab_datenum(d.year, d.month, d.day, d.hour, d.minute, d.second)
    # datetime.utcfromtimestamp(Date_Num[i])

# Time series of POC flux in the eddy core
width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.4, width, height - 0.15])
plt.plot(x_filtered, Flux_eddy_core)
plt.ylabel('Flux(mgC/$m^2$/d)')
plt.ylim(0, ax.get_ylim()[1])
plt.xlim(x_filtered[0], x_filtered[-1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(x_filtered.min(), x_filtered.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# I add the panel label
ax.text(-0.05, 1.05, 'f', transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# I save
plt.savefig('../Plots/Fig_Main_v07/Supplementary/TimeSeries_EddyCore/TimeSeries_EddyCore_Flux_v07.pdf' , dpi=200)
plt.close()

#######################################################################
# I save some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=68;print(matlab_datevec(x_filtered_datenum[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0413to0803_eddy_core_down'
arg_value=np.mean(Flux_depthf[0:i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Flux_0413to0803_eddy_core_up'
arg_value=np.mean(Flux_depth0[0:i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=73;print(matlab_datevec(x_filtered_datenum[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0812_eddy_core_up'
arg_value=np.mean(Flux_depth0[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i1=67;print(matlab_datevec(x_filtered_datenum[i1]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
i2=86;print(matlab_datevec(x_filtered_datenum[i2]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_August_eddy_core_up'
arg_value=np.mean(Flux_depth0[i1:i2])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=44;print(matlab_datevec(mld_datenum0[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0810_ML'
arg_value=np.mean(Flux_MLD[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=39;print(matlab_datevec(mld_datenum0[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'MajorFlux_ML_start'
write_latex_data(filename,argument,'26 July')
i=53;print(matlab_datevec(mld_datenum0[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'MajorFlux_ML_end'
write_latex_data(filename,argument,'8 September')
argument = 'FluxMiP_eddy_core_up_percentage'
arg_value=np.mean(Flux_MiP_depth0/(Flux_MiP_depth0+Flux_MaP_depth0)*100)
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'FluxMiP_eddy_core_down_percentage'
arg_value=np.mean(Flux_MiP_depthf/(Flux_MiP_depthf+Flux_MaP_depthf)*100)
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'Flux_0413_eddy_core'
arg_value=np.mean(Flux_eddy_core[0])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=61;print(matlab_datevec(x_filtered_datenum[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0723_eddy_core'
arg_value=np.mean(Flux_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=84;print(matlab_datevec(x_filtered_datenum[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0830_eddy_core'
arg_value=np.mean(Flux_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)


#endregion


#######################################################################################################################
######### Fig. 03b
########################################################################################################################
# region Fig. 03b
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
#######
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to finish the time series in the plot
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
delta_depth=15                  # around of the depth which I consider when extracting the flux

#######################################################################
# I load the ecopart data
#######################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data_ecopart=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data_ecopart.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data
lon=np.array(data_ecopart['Longitude'][sel_filename])
lat=np.array(data_ecopart['Latitude'][sel_filename])
Date_Time=np.array(data_ecopart['Date_Time'][sel_filename])
Flux=np.array(data_ecopart['Flux_mgC_m2'][sel_filename])
depth=np.array(data_ecopart['Depth [m]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[0:sel_insideEddy.size];list_dates=list_dates[sel_insideEddy]
# list_dates=list_dates[sel_insideEddy[0:list_dates.size]]

#######################################################################
# I filter and interpolate the flux
#######################################################################

day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

parameter=Flux.copy()
parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
# I filter the flux prophiles
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
    z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        parameter_filtered = np.concatenate((parameter_filtered, z))
        Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
        depth_filtered = np.concatenate((depth_filtered, y2))

x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
# I interpolate
parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

depthf = 600  # final depth
sel_depthf_600 = (np.abs(y_filtered) > depthf - delta_depth) & (np.abs(y_filtered) <= depthf + delta_depth)
Flux_interp_600 = np.mean(parameter_interp[sel_depthf_600, :], axis=0)

depthf = 200  # final depth
sel_depthf_200 = (np.abs(y_filtered) > depthf - delta_depth) & (np.abs(y_filtered) <= depthf + delta_depth)
Flux_interp_200 = np.mean(parameter_interp[sel_depthf_200, :], axis=0)



########################################################################################################################
######### Plot Fig. 03b
########################################################################################################################
max_Flux=32

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
parameter_plot=parameter_interp.copy()
parameter_plot[parameter_plot<0]=0
parameter_plot[parameter_plot>max_Flux]=max_Flux
ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
plt.gca().invert_yaxis()
plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
# I draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the panel label
ax.text(-0.05, 1.05, 'b', transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig03b_v07.pdf' ,dpi=200)
plt.close()

#######################################################################
# I save some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=68;print(matlab_datevec(x_filtered_datenum[i]).astype(int)) #Day in which flux at upper eddy core becomes larger than the one at lower eddy core
argument = 'Flux_0413to0803_eddy_core_low'
arg_value=np.mean(Flux_depthf[0:i])
write_latex_data(filename,argument,'%0.2f' % arg_value)

# endregion

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 04: Carbon budget. Fig. Suppl.: Carbon budget extended. Fig. Suppl.: Carbon budget , second TW. Fig. Suppl.: Carbon budget extended, second TW.
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Fig. 04
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
########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from paruvpy import bbpPOC_to_bbp_abundance
from paruvpy import calculate_remin_func_abun
from paruvpy import ESD_limits_to_ESD_middle
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
import seawater as sw
import gsw
from lin_fit import lin_fit

#######################################################################
# I define the function to extract the POC content at the beginning and end of the time series, used in the function below
#######################################################################
def POC_day0dayf(t,POC):
    ndays=t.max()
    (interpol, slpe_ci, intrcp_ci, _, _) = lin_fit(t, POC)
    POC_day0=interpol.intercept / ndays# * layer_thickness
    POC_dayf=(interpol.slope*ndays+interpol.intercept) / ndays# * layer_thickness
    slpe_std=abs((slpe_ci - interpol.slope)[0]);intrcp_std=abs((intrcp_ci - interpol.intercept)[0])
    POC_day0_std = intrcp_std / ndays# * layer_thickness
    POC_dayf_std = np.sqrt(ndays**2*slpe_std**2+intrcp_std**2) / ndays# * layer_thickness
    return POC_day0,POC_dayf,POC_day0_std,POC_dayf_std

#######################################################################
# I define the function for the carbon budget calculation
#######################################################################
# region carbon_budget_calculation(dens0,densf,day0,dayf):
def carbon_budget_calculation(dens0,densf,day0,dayf):
    ########################################################################################################################
    # Starting parameters
    ########################################################################################################################
    # dayf = day0+timedelta(days=ndays) # final date for the carbon budget calculation
    ndays = (dayf - day0).days  # number of days
    delta_dens_flux = 0.025     # around of the density which I consider when extracting the flux
    Oxy2C = 0.89                # to convert from mol of oxygen to mol of carbon
    Oxy2C_std = Oxy2C * 0.4
    mol2gC = 12.0107            # to convert from mol of carbon to grams of carbon
    day0_float = calendar.timegm(day0.timetuple())
    dayf_float = calendar.timegm(dayf.timetuple())
    day0_datenum = matlab_datenum(day0.year,day0.month,day0.day,day0.hour,day0.minute,day0.second)
    dayf_datenum = matlab_datenum(dayf.year,dayf.month,dayf.day,dayf.hour,dayf.minute,dayf.second)

    ########################################################################################################################
    # I load and process data
    ########################################################################################################################

    #I load the file with the flux and POC
    filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
    data=pd.read_csv(filename_ecopart, sep='\t', header=0)
    RAWfilename=data.RAWfilename

    #I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
    ct=0
    sel_filename = [True for i in range(RAWfilename.size)]
    for a in RAWfilename:
        if a.split('-')[-1].split('_')[0] == 'ASC':
            sel_filename[ct]=True
        else:
            sel_filename[ct] = False
        ct+=1

    # I extract the data
    Date_Time=np.array(data['Date_Time'][sel_filename])
    depth=np.array(data['Depth [m]'][sel_filename])
    dens=np.array(data['Potential density [kg/m3]'][sel_filename])
    Flux=np.array(data['Flux_mgC_m2'][sel_filename])
    # Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
    Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
    # Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
    MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
    MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
    MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num=np.r_[0:Flux.size]
    for i in Date_Num:
        date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    # I select the data only in the prescribed period
    list_dates=np.sort(np.unique(Date_Num))
    list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
    n_profiles=list_dates.size

    # I load the bbp data and I select only those in the prescribed period
    storedir = '%s/GIT/AC_Agulhas_eddy_2021/Data' % home
    a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
    data_an18 = pickle.load(a_file)
    bbp_POC = data_an18['bbp_POC']
    bbp_POC_Koestner = data_an18['bbp_POC_Koestner']
    Date_Num_bbp = data_an18['Date_Num_bbp']
    Date_Vec_bbp = data_an18['Date_Vec_bbp']
    depth_bbp = data_an18['depth_bbp']
    dens_bbp = data_an18['dens_bbp']
    temp_bbp = data_an18['temperature']
    a_file.close()

    sel_dates = (Date_Num_bbp>=day0_datenum)&(Date_Num_bbp<=dayf_datenum)
    Date_Num_bbp = Date_Num_bbp[sel_dates]
    Date_Vec_bbp = Date_Vec_bbp[sel_dates,:]
    depth_bbp = depth_bbp[sel_dates,:]
    dens_bbp = dens_bbp[sel_dates,:]
    temp_bbp = temp_bbp[sel_dates,:]
    bbp_POC = bbp_POC[sel_dates, :]
    bbp_POC_Koestner = bbp_POC_Koestner[sel_dates, :]
    # I calculate the bbp abundance
    sel=(bbp_POC==99999)&(bbp_POC_Koestner==99999)
    bbp_abun_tmp=bbpPOC_to_bbp_abundance(bbp_POC[~sel])
    bbp_abun=bbp_POC.copy()*0+99999;bbp_abun[~sel]=bbp_abun_tmp
    bbp_abun_Koestner_tmp=bbpPOC_to_bbp_abundance(bbp_POC_Koestner[~sel])
    bbp_abun_Koestner=bbp_POC.copy()*0+99999;bbp_abun_Koestner[~sel]=bbp_abun_Koestner_tmp
    # I calculate the bbp PARR in nmol_l_h
    esd_bbp = ESD_limits_to_ESD_middle(0.001, 0.025)
    bbp_PARR_tmp = calculate_remin_func_abun(bbp_abun[~sel], esd_bbp, temp_bbp[~sel], kRemPoc=0.013)
    bbp_PARR=bbp_POC.copy()*0+99999;bbp_PARR[~sel]=bbp_PARR_tmp
    bbp_PARR_Koestner_tmp = calculate_remin_func_abun(bbp_abun_Koestner[~sel], esd_bbp, temp_bbp[~sel], kRemPoc=0.013)
    bbp_PARR_Koestner=bbp_POC.copy()*0+99999;bbp_PARR_Koestner[~sel]=bbp_PARR_Koestner_tmp
    del bbp_abun_tmp,bbp_abun_Koestner_tmp,bbp_abun,bbp_abun_Koestner,bbp_PARR_tmp,bbp_PARR_Koestner_tmp
    # I convert the bpp PARR from nmol_l_h to micromol/kg/day
    bbp_PARR[~sel] = bbp_PARR[~sel] / 1000 * 24 / (dens_bbp[~sel] / 1000)
    bbp_PARR_Koestner[~sel] = bbp_PARR_Koestner[~sel] / 1000 * 24 / (dens_bbp[~sel] / 1000)

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_bbp_calendar = Date_Num_bbp.copy()
    for i in range(0,Date_Num_bbp_calendar.size):
        date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                                          Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
        Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
        # datetime.utcfromtimestamp(Date_Num[i])

    ########################################################################################################################
    # Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because
    # some profiles only reach 400 m; (iii) I extract the mean MiP+MaP+bbp values between dens0 and densf and between day0 and
    # dayf (I obtain a time series); (iv) I calculate the bbp PARR between dens0 and densf and between day0 and dayf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
    MiP_extended_filtered=np.array([]);dens_MiP_extended_filtered=np.array([]);Date_Num_MiP_extended_filtered=np.array([])
    MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
    bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])
    bbp_Koestner_filtered=np.array([]);dens_bbp_Koestner_filtered=np.array([]);Date_Num_bbp_Koestner_filtered=np.array([])
    bbp_PARR_filtered=np.array([]);dens_bbp_PARR_filtered=np.array([]);Date_Num_bbp_PARR_filtered=np.array([])
    bbp_PARR_Koestner_filtered=np.array([]);dens_bbp_PARR_Koestner_filtered=np.array([]);Date_Num_bbp_PARR_Koestner_filtered=np.array([])

    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=MiP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_filtered = np.concatenate((MiP_filtered, z))
            Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x2))
            dens_MiP_filtered = np.concatenate((dens_MiP_filtered, y2))
        z=MiP_POC_extended[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_extended_filtered = np.concatenate((MiP_extended_filtered, z))
            Date_Num_MiP_extended_filtered = np.concatenate((Date_Num_MiP_extended_filtered, x2))
            dens_MiP_extended_filtered = np.concatenate((dens_MiP_extended_filtered, y2))
        z=MaP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MaP_filtered = np.concatenate((MaP_filtered, z))
            Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x2))
            dens_MaP_filtered = np.concatenate((dens_MaP_filtered, y2))

    i=0
    for i in range(0, bbp_POC.shape[0]):
        #Cetinic
        z=bbp_POC[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>100] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_filtered = np.concatenate((bbp_filtered, z))
            Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_filtered = np.concatenate((dens_bbp_filtered, y2))
        #Cetinic PARR
        z=bbp_PARR[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>100] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_PARR_filtered = np.concatenate((bbp_PARR_filtered, z))
            Date_Num_bbp_PARR_filtered = np.concatenate((Date_Num_bbp_PARR_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_PARR_filtered = np.concatenate((dens_bbp_PARR_filtered, y2))
        #Koestner
        z=bbp_POC_Koestner[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>400] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_Koestner_filtered = np.concatenate((bbp_Koestner_filtered, z))
            Date_Num_bbp_Koestner_filtered = np.concatenate((Date_Num_bbp_Koestner_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_Koestner_filtered = np.concatenate((dens_bbp_Koestner_filtered, y2))
        #Koestner PARR
        z=bbp_PARR_Koestner[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>400] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_PARR_Koestner_filtered = np.concatenate((bbp_PARR_Koestner_filtered, z))
            Date_Num_bbp_PARR_Koestner_filtered = np.concatenate((Date_Num_bbp_PARR_Koestner_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_PARR_Koestner_filtered = np.concatenate((dens_bbp_PARR_Koestner_filtered, y2))

    # I define the x and y arrays for the MiP+MaP+bbp interpolation
    x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), ndays)
    y_filtered = np.linspace(dens_bbp_filtered.min(), dens_MaP_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    MiP_interp = griddata((Date_Num_MiP_filtered, dens_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MiP_extended_interp = griddata((Date_Num_MiP_extended_filtered, dens_MiP_extended_filtered), MiP_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MaP_interp = griddata((Date_Num_MaP_filtered, dens_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_interp = griddata((Date_Num_bbp_filtered, dens_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_Koestner_interp = griddata((Date_Num_bbp_Koestner_filtered, dens_bbp_Koestner_filtered), bbp_Koestner_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_PARR_interp = griddata((Date_Num_bbp_PARR_filtered, dens_bbp_PARR_filtered), bbp_PARR_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_PARR_Koestner_interp = griddata((Date_Num_bbp_PARR_Koestner_filtered, dens_bbp_PARR_Koestner_filtered), bbp_PARR_Koestner_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, I calculate the mean MiP+MaP+bbp (and std) between dens0 and densf between day0 and dayf
    sel_dens0_densf = (np.abs(y_filtered) >= dens0) & (np.abs(y_filtered) < densf)
    MiP_POC_dens0_densf=np.mean(MiP_interp[sel_dens0_densf,:],0)
    MiP_POC_extended_dens0_densf=np.mean(MiP_extended_interp[sel_dens0_densf,:],0)
    MaP_POC_dens0_densf=np.mean(MaP_interp[sel_dens0_densf,:],0)
    bbp_POC_dens0_densf=np.mean(bbp_interp[sel_dens0_densf,:],0)
    bbp_POC_Koestner_dens0_densf=np.mean(bbp_Koestner_interp[sel_dens0_densf,:],0)
    bbp_PARR_dens0_densf=np.mean(bbp_PARR_interp[sel_dens0_densf,:])
    bbp_PARR_Koestner_dens0_densf=np.mean(bbp_PARR_Koestner_interp[sel_dens0_densf,:])

    MiP_POC_dens0_densf_std = np.std(MiP_interp[sel_dens0_densf, :], 0)
    MiP_POC_extended_dens0_densf_std = np.std(MiP_extended_interp[sel_dens0_densf, :], 0)
    MaP_POC_dens0_densf_std = np.std(MaP_interp[sel_dens0_densf, :], 0)
    bbp_POC_dens0_densf_std = np.std(bbp_interp[sel_dens0_densf, :], 0)
    bbp_POC_Koestner_dens0_densf_std = np.std(bbp_Koestner_interp[sel_dens0_densf, :], 0)
    bbp_PARR_dens0_densf_std = np.std(bbp_PARR_interp[sel_dens0_densf, :])
    bbp_PARR_Koestner_dens0_densf_std = np.std(bbp_PARR_Koestner_interp[sel_dens0_densf, :])

    Integrated_POC_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    Integrated_POC_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    Integrated_POC_noBBP_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf
    Integrated_POC_noBBP_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf
    Integrated_POC_noBBP_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2  )
    Integrated_POC_noBBP_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 )
    Integrated_POC_Koestner_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf + bbp_POC_Koestner_dens0_densf
    Integrated_POC_Koestner_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf + bbp_POC_Koestner_dens0_densf
    Integrated_POC_Koestner_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_Koestner_dens0_densf_std**2 )
    Integrated_POC_Koestner_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_Koestner_dens0_densf_std**2 )

    ########################################################################################################################
    # Here I extract the flux values at dens0 and densf. To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid in time and density. This step is necessary to have the flux at 600 m, because some
    # profiles only reach 400 m; (iii) I extract the flux values at dens0 and densf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    Flux_filtered=np.array([]);dens_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
    Flux_extended_filtered=np.array([]);dens_Flux_extended_filtered=np.array([]);Date_Num_Flux_extended_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=Flux[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_filtered = np.concatenate((Flux_filtered, z))
            Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
            dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y2))
        z=Flux_extended[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered = np.concatenate((Flux_extended_filtered, z))
            Date_Num_Flux_extended_filtered = np.concatenate((Date_Num_Flux_extended_filtered, x2))
            dens_Flux_extended_filtered = np.concatenate((dens_Flux_extended_filtered, y2))

    # I define the x and y arrays for the Flux interpolation
    x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
    y_filtered = np.linspace(dens_Flux_filtered.min(), dens_Flux_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    Flux_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_interp = griddata((Date_Num_Flux_extended_filtered, dens_Flux_extended_filtered), Flux_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, flux extraction at dens0 and densf

    sel_layer = (np.abs(y_filtered) >= dens0-delta_dens_flux) & (np.abs(y_filtered) < dens0+delta_dens_flux)
    Flux_dens0 = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_extended_dens0 = np.mean(Flux_extended_interp[sel_layer, :], axis=0)

    sel_layer = (np.abs(y_filtered) >= densf - delta_dens_flux) & (np.abs(y_filtered) < densf + delta_dens_flux)
    Flux_densf = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_extended_densf = np.mean(Flux_extended_interp[sel_layer,:],axis=0)


    ########################################################################################################################
    # Here I calculate the carbon consumption rate due to (i) oxygen consumption and (ii) PARR
    ########################################################################################################################

    ############### I load Coriolis data with the oxygen information
    filename='6903095_Sprof_all.nc'
    ds = nc.Dataset('%s/%s' % (storedir,filename))

    lon=np.array(ds.variables['LONGITUDE'])
    lat=np.array(ds.variables['LATITUDE'])
    Date_Num=np.array(ds.variables['JULD'])
    temp=np.array(ds.variables['TEMP_ADJUSTED'])
    pres=np.array(ds.variables['PRES_ADJUSTED'])
    psal=np.array(ds.variables['PSAL_ADJUSTED'])
    doxy=np.array(ds.variables['DOXY_ADJUSTED'])

    #If adjusted values are not available yet, I take the non adjusted ones
    if np.sum(temp==99999)==temp.size:
        print('Taking non adjusted temperature')
        temp = np.array(ds.variables['TEMP'])
        temp_qc = np.array(ds.variables['TEMP_QC'])
    if np.sum(pres==99999)==pres.size:
        print('Taking non adjusted pressure')
        pres = np.array(ds.variables['PRES'])
        pres_qc = np.array(ds.variables['PRES_QC'])
    if np.sum(psal==99999)==psal.size:
        print('Taking non adjusted salinity')
        psal = np.array(ds.variables['PSAL'])
        psal_qc = np.array(ds.variables['PSAL_QC'])
    if np.sum(doxy==99999)==doxy.size:
        print('Taking non adjusted oxygen')
        doxy = np.array(ds.variables['DOXY'])
        doxy_qc = np.array(ds.variables['DOXY_QC'])

    #I tranform the pressure to depth
    mask_depth=pres!=99999 #I select only valid values
    lat_tmp=np.tile(lat,[pres.shape[1],1]).T
    lat_tmp=lat_tmp[mask_depth]
    pres_tmp=pres[mask_depth]
    depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
    depth=np.ones(temp.shape)*99999
    depth[mask_depth]=depth_tmp

    #I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
    #salinity and temperature first
    mask_dens=np.logical_and(pres!=99999,temp!=99999,psal!=99999) # I exclude the points with value = 99999
    lat_tmp=np.tile(lat,[pres.shape[1],1]).T
    lon_tmp=np.tile(lon,[pres.shape[1],1]).T
    lat_tmp=lat_tmp[mask_dens]
    lon_tmp=lon_tmp[mask_dens]
    pres_tmp=pres[mask_dens]
    psal_tmp=psal[mask_dens]
    temp_tmp=temp[mask_dens]
    abs_psal_tmp=gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp) # I compute absolute salinity
    cons_tmp=gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)          # I compute conservative temperature
    dens_tmp=gsw.density.sigma0(abs_psal_tmp, cons_tmp)
    dens=np.ones(temp.shape)*99999
    dens[mask_dens]=dens_tmp+1000

    ############### I load the PARR correlated data (density, depth, and date) from Ecopart
    Date_Time_PARR=np.array(data['Date_Time'][sel_filename])
    depth_PARR=np.array(data['Depth [m]'][sel_filename])
    dens_PARR=np.array(data['Potential density [kg/m3]'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_PARR=np.r_[0:Date_Time_PARR.size]
    for i in range(0,Date_Time_PARR.size):
        date_time_obj = datetime.datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    list_dates_PARR=np.unique(Date_Num_PARR)

    #############################################
    ############### Loop on the different isopycnal values chosen for the study of the oxygen profile
    Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+ndays]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))

    # Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration between
    # dens0 and densf
    reference_isopycnal=(dens0+densf)*0.5
    reference_isopycnal_down=dens0
    reference_isopycnal_up=densf
    doxy_isopycnal=np.array([]);depth_isopycnal_tmp=np.array([]);depth_isopycnal_down_tmp=np.array([]);depth_isopycnal_up_tmp=np.array([])
    Date_Num_isopycnal=np.array([])
    i=0
    for i in range(0,doxy.shape[0]):
        #Here, for the i-th profile, I select the oxygen, density and depth profiles of Coriolis data, excluding the nan values
        sel = (doxy[i, :] != 99999) & (dens[i, :] != 99999)
        z=doxy[i,sel];y=dens[i,sel];d=depth[i,sel]

        # Here I proceed only if the date is inside the Date_Num_limit fixed
        if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

            # Here I extract the oxygen along the isopycnal
            sel_layer = (y >= reference_isopycnal_down) & (y < reference_isopycnal_up)
            if np.sum(sel_layer) > 0:  # If sel_layer has some True values, then I take as the doxy of this isopycnal the mean of the doxy values in correspondence with these layers
                doxy_tmp = np.mean(z[sel_layer])
                depth_isopycnal_tmp2 = d[sel_layer]
                doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, np.mean(depth_isopycnal_tmp2) )
                depth_isopycnal_down_tmp = np.append(depth_isopycnal_down_tmp, depth_isopycnal_tmp2[0] )
                depth_isopycnal_up_tmp = np.append(depth_isopycnal_up_tmp, depth_isopycnal_tmp2[-1] )
            else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the oxygen concentration at that depth
                depth_isopycnal_tmp2 = np.array([])
                for iy in range(0, y.size - 1):
                    if y[iy] <= reference_isopycnal < y[iy + 1]:
                        dist = (reference_isopycnal - y[iy]) / (y[iy + 1] - y[iy])
                        doxy_tmp = z[iy] + (z[iy + 1] - z[iy]) * dist
                        d_tmp = d[iy] + (d[iy + 1] - d[iy]) * dist
                        doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                        depth_isopycnal_tmp2 = np.append(depth_isopycnal_tmp2, d_tmp)
                        Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                if depth_isopycnal_tmp2.size > 0:
                    depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, np.mean(depth_isopycnal_tmp2))
                    depth_isopycnal_down_tmp = np.append(depth_isopycnal_down_tmp, depth_isopycnal_tmp2[0])
                    depth_isopycnal_up_tmp = np.append(depth_isopycnal_up_tmp, depth_isopycnal_tmp2[-1])

    # If I have at least three points I interpolate them, so that the slope gives me the respiration rate
    # (in micromol/kg per day)
    depth_isopycnal = np.mean(depth_isopycnal_tmp)
    depth_isopycnal_down = np.mean(depth_isopycnal_down_tmp)
    depth_isopycnal_up = np.mean(depth_isopycnal_up_tmp)
    layer_thickness=np.mean(depth_isopycnal_up_tmp - depth_isopycnal_down_tmp)
    if Date_Num_isopycnal.size>2:
        (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)
        # fig = plt.figure(1, figsize=(12, 8))
        # plot1 = plt.scatter(Date_Num_isopycnal,doxy_isopycnal)
        slope_doxy = interpol.slope
        slope_ci_doxy = np.reshape(slpe_ci.copy(),(1,2))

        O2_resp_mgO_m3_d=-slope_doxy.copy()*reference_isopycnal*mol2gC/1000#*ndays*layer_thickness
        O2_resp_mgO_m3_d_ci=-slope_ci_doxy.copy()*( np.tile(reference_isopycnal,(2,1)).T) *mol2gC/1000#*ndays*layer_thickness
        O2_resp_mgC_m3_d = O2_resp_mgO_m3_d*Oxy2C
        #Since I have an uncertainty both on O2_resp_mgO_m3_d and on Oxy2C, I propagate the error to obtain the uncertainty on O2_resp_mgC_m3_d
        O2_resp_mgO_m3_d_ci=abs(np.diff(O2_resp_mgO_m3_d_ci)[0][0]/2)
        O2_resp_mgC_m3_d_err = np.sqrt( Oxy2C**2*O2_resp_mgO_m3_d_ci**2 + Oxy2C_std**2*O2_resp_mgO_m3_d**2 )
        O2_resp_mgC_m3_d_ci=np.ones((1,2));O2_resp_mgC_m3_d_ci[0,0]=O2_resp_mgC_m3_d-O2_resp_mgC_m3_d_err;O2_resp_mgC_m3_d_ci[0,1]=O2_resp_mgC_m3_d+O2_resp_mgC_m3_d_err

    else:
        O2_resp_mgC_m3_d=np.nan
        O2_resp_mgC_m3_d_ci=np.nan*np.ones((1,2))

    #############################################
    ############### Loop on different respiration types used to estimate PARR
    #List of the different Respiration types present in data
    list_Respi_types = [match for match in data.columns if "Respi" in match]
    nRespi= len(list_Respi_types)  #number of respiration types

    POC_resp_mgC_m3_d_list=np.zeros((nRespi,))
    POC_resp_mgC_m3_d_std_list = np.zeros((nRespi,))

    iRespi=0
    for iRespi in range(0,nRespi):
        PARR_nmol_l_h = np.array(data[list_Respi_types[iRespi]][sel_filename])
        # I convert the PARR measured in micromol/kg/day
        PARR_micromol_kg_day = PARR_nmol_l_h.copy() / 1000 * 24 / (dens_PARR / 1000)

        #############################################
        ############### Loop on the different PARR profiles

        #Here, for each profile included between two dates (Date_Num_limit), I compute PARR concentration between dens0
        # and densf
        PARR_isopycnal=np.array([]);#depth_PARR_tmp=np.array([]);dens_PARR_tmp=np.array([])
        i=0
        for i in range(0,doxy.shape[0]):
            # Here, for the i-th profile, I select the PARR, density and depth profiles of Ecopart data, excluding the nan values. I do the same for the bbp
            sel_PARR=Date_Num_PARR==list_dates_PARR[i]
            z_PARR=PARR_micromol_kg_day[sel_PARR];y_PARR=dens_PARR[sel_PARR];d_PARR=depth_PARR[sel_PARR]
            sel_PARR = (~np.isnan(z_PARR)) & (~np.isnan(y_PARR))
            z_PARR = z_PARR[sel_PARR];y_PARR = y_PARR[sel_PARR];d_PARR = d_PARR[sel_PARR]

            # Here I proceed only if the date is inside the Date_Num_limit fixed
            if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

                # Here I extract the PARR along the isopycnal
                sel_layer_PARR = (y_PARR >= reference_isopycnal_down) & (y_PARR < reference_isopycnal_up)
                if np.sum(sel_layer_PARR) > 0:  # If sel_layer_PARR has some True values, then I take as the PARR of this isopycnal the mean of the PARR values in correspondence with these layers
                    PARR_tmp = np.mean(z_PARR[sel_layer_PARR])
                    PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the PARR at that depth
                    for iy in range(0, y_PARR.size - 1):
                        if y_PARR[iy] <= reference_isopycnal < y_PARR[iy + 1]:
                            dist = (reference_isopycnal - y_PARR[iy]) / (y_PARR[iy + 1] - y_PARR[iy])
                            PARR_tmp = z_PARR[iy] + (z_PARR[iy + 1] - z_PARR[iy]) * dist
                            d_tmp = d_PARR[iy] + (d_PARR[iy + 1] - d_PARR[iy]) * dist
                            PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)

        PARR_isopycnal_std = np.std(PARR_isopycnal)
        PARR_isopycnal = np.mean(PARR_isopycnal)

        # I convert the PARR and the oxygen respiration rates (in micromolO2/kg/d) to the total amount of carbon consumption
        # between depth0 and depthf, and between day0 and dayf (in mgC/m3/d)
        # *Oxy2C -> to micromolC/kg/d
        # *mol2gC -> to microgC/kg/d
        # /1000 -> to mgC/kg/d
        # *density -> to mgC/m3/d
        # *layer_thickness*ndays -> to mgC/m2

        # I calculate the PARR
        POC_resp_mgC_m3_d_list[iRespi] = PARR_isopycnal.copy()*reference_isopycnal*Oxy2C*mol2gC/1000 # *ndays * layer_thickness
        # I calculate the error on the PARR (POC_resp) by taking into account the error on Oxy2C as well
        POC_resp_mgC_m3_d_std_list[iRespi] = reference_isopycnal*mol2gC/1000*np.sqrt(PARR_isopycnal.copy()**2*Oxy2C_std**2+PARR_isopycnal_std.copy()**2*Oxy2C**2) # *ndays * layer_thickness

    # I convert bbp to mgC_m3_d: first, I calculate the error otherwise the formula is wrong
    # I calculate the error on the bbp_PARR by taking into account the error on Oxy2C as well
    bbp_PARR_dens0_densf_std = reference_isopycnal * mol2gC / 1000*np.sqrt(bbp_PARR_dens0_densf_std.copy()**2*Oxy2C**2 +bbp_PARR_dens0_densf.copy()**2*Oxy2C_std**2 ) # *ndays * layer_thickness
    bbp_PARR_Koestner_dens0_densf_std = reference_isopycnal * mol2gC / 1000*np.sqrt(bbp_PARR_Koestner_dens0_densf_std.copy()**2*Oxy2C**2 +bbp_PARR_Koestner_dens0_densf.copy()**2*Oxy2C_std**2 )   # *ndays * layer_thickness
    # I convert bbp to mgC_m3_d
    bbp_PARR_dens0_densf = bbp_PARR_dens0_densf.copy() * reference_isopycnal * Oxy2C * mol2gC / 1000  # *ndays * layer_thickness
    bbp_PARR_Koestner_dens0_densf = bbp_PARR_Koestner_dens0_densf.copy() * reference_isopycnal * Oxy2C * mol2gC / 1000  # *ndays * layer_thickness

    ########################################################################################################################
    # Here I calculate the carbon budget for depth0—depthf layer
    ########################################################################################################################
    # Date_Num_Flux = x_filtered
    # depth_POC_resp = list_depth_PARR
    # depth_02_resp = depth_isopycnal

    ############### I calculate the integrated POC (MiP+MaP+bbp), between depth0 and depthf, for day0 and dayf. I transform it to mgC/m2

    # t=np.r_[0:ndays]
    # (Integrated_POC_day0_mgC_m3_d,Integrated_POC_dayf_mgC_m3_d,Integrated_POC_day0_mgC_m3_d_std,Integrated_POC_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_mgC_m3)
    # (Integrated_POC_noBBP_day0_mgC_m3_d,Integrated_POC_noBBP_dayf_mgC_m3_d,Integrated_POC_noBBP_day0_mgC_m3_d_std,Integrated_POC_noBBP_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_noBBP_mgC_m3)
    # (Integrated_POC_Koestner_day0_mgC_m3_d,Integrated_POC_Koestner_dayf_mgC_m3_d,Integrated_POC_Koestner_day0_mgC_m3_d_std,Integrated_POC_Koestner_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_Koestner_mgC_m3)

    # I calculate the integrated POC between day and dayf, based on using the average value between day0 and day0+3days and between dayf-3days and dayf, respectively
    nd = 3
    Integrated_POC_day0_mgC_m3_d = np.nanmean(Integrated_POC_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_dayf_mgC_m3_d = np.nanmean(Integrated_POC_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_day0_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_noBBP_dayf_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_noBBP_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_day0_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_Koestner_dayf_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_Koestner_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Delta_Integrated_POC = Integrated_POC_dayf_mgC_m3_d - Integrated_POC_day0_mgC_m3_d
    Delta_Integrated_POC_std = np.sqrt( Integrated_POC_dayf_mgC_m3_d_std**2 + Integrated_POC_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_noBBP = Integrated_POC_noBBP_dayf_mgC_m3_d - Integrated_POC_noBBP_day0_mgC_m3_d
    Delta_Integrated_POC_noBBP_std = np.sqrt( Integrated_POC_noBBP_dayf_mgC_m3_d_std**2 + Integrated_POC_noBBP_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_Koestner = Integrated_POC_Koestner_dayf_mgC_m3_d - Integrated_POC_Koestner_day0_mgC_m3_d
    Delta_Integrated_POC_Koestner_std = np.sqrt( Integrated_POC_Koestner_dayf_mgC_m3_d_std**2 + Integrated_POC_Koestner_day0_mgC_m3_d_std**2 )

    # (Integrated_POC_extended_day0_mgC_m3_d,Integrated_POC_extended_dayf_mgC_m3_d,Integrated_POC_extended_day0_mgC_m3_d_std,Integrated_POC_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_extended_mgC_m3)
    # (Integrated_POC_noBBP_extended_day0_mgC_m3_d,Integrated_POC_noBBP_extended_dayf_mgC_m3_d,Integrated_POC_noBBP_extended_day0_mgC_m3_d_std,Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_noBBP_extended_mgC_m3)
    # (Integrated_POC_Koestner_extended_day0_mgC_m3_d,Integrated_POC_Koestner_extended_dayf_mgC_m3_d,Integrated_POC_Koestner_extended_day0_mgC_m3_d_std,Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_Koestner_extended_mgC_m3)

    # Old way of calculate the integrated POC between day and dayf, based on using the exact value on day0 and dayf, respectively
    Integrated_POC_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Delta_Integrated_POC_extended = Integrated_POC_extended_dayf_mgC_m3_d - Integrated_POC_extended_day0_mgC_m3_d
    Delta_Integrated_POC_extended_std = np.sqrt( Integrated_POC_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_extended_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_noBBP_extended = Integrated_POC_noBBP_extended_dayf_mgC_m3_d - Integrated_POC_noBBP_extended_day0_mgC_m3_d
    Delta_Integrated_POC_noBBP_extended_std = np.sqrt( Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_noBBP_extended_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_Koestner_extended = Integrated_POC_Koestner_extended_dayf_mgC_m3_d - Integrated_POC_Koestner_extended_day0_mgC_m3_d
    Delta_Integrated_POC_Koestner_extended_std = np.sqrt( Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_Koestner_extended_day0_mgC_m3_d_std**2 )

    ############### I calculate the amount of POC entering from depht0 and exiting from dayf between day0 and dayf (in mgC/m3/day)

    # I extract the index of Flux_dens0/Flux_densf which correspond to day0 (and dayf)
    Flux_dens0_mgC_m3_d = np.mean(Flux_dens0) / layer_thickness # * ndays
    Flux_dens0_mgC_m3_d_std = np.std(Flux_dens0) / layer_thickness # * ndays
    Flux_densf_mgC_m3_d = np.mean(Flux_densf) / layer_thickness # * ndays
    Flux_densf_mgC_m3_d_std = np.std(Flux_densf) / layer_thickness # * ndays

    Flux_extended_dens0_mgC_m3_d = np.mean(Flux_extended_dens0) / layer_thickness # * ndays
    Flux_extended_dens0_mgC_m3_d_std = np.std(Flux_extended_dens0) / layer_thickness # * ndays
    Flux_extended_densf_mgC_m3_d = np.mean(Flux_extended_densf) / layer_thickness # * ndays
    Flux_extended_densf_mgC_m3_d_std = np.std(Flux_extended_densf) / layer_thickness # * ndays

    Delta_flux = Flux_dens0_mgC_m3_d - Flux_densf_mgC_m3_d
    Delta_flux_extended = Flux_extended_dens0_mgC_m3_d - Flux_extended_densf_mgC_m3_d

    Delta_flux_std = np.sqrt( Flux_dens0_mgC_m3_d_std**2 + Flux_densf_mgC_m3_d_std**2 )
    Delta_flux_extended_std = np.sqrt( Flux_extended_dens0_mgC_m3_d_std**2 + Flux_extended_densf_mgC_m3_d_std**2 )

    Theoretical_Budget = Delta_flux - Delta_Integrated_POC
    Theoretical_Budget_extended = Delta_flux_extended - Delta_Integrated_POC_extended
    Theoretical_Budget_noBBP = Delta_flux - Delta_Integrated_POC_noBBP
    Theoretical_Budget_noBBP_extended = Delta_flux_extended - Delta_Integrated_POC_noBBP_extended
    Theoretical_Budget_Koestner = Delta_flux - Delta_Integrated_POC_Koestner
    Theoretical_Budget_Koestner_extended = Delta_flux_extended - Delta_Integrated_POC_Koestner_extended

    Theoretical_Budget_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_extended_std**2 )
    Theoretical_Budget_noBBP_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_noBBP_std**2 )
    Theoretical_Budget_noBBP_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_noBBP_extended_std**2 )
    Theoretical_Budget_Koestner_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_Koestner_std**2 )
    Theoretical_Budget_Koestner_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_Koestner_extended_std**2 )


    ############### I return the data
    return Theoretical_Budget,Theoretical_Budget_std, \
           Theoretical_Budget_extended,Theoretical_Budget_extended_std, \
           Theoretical_Budget_noBBP,Theoretical_Budget_noBBP_std, \
           Theoretical_Budget_noBBP_extended,Theoretical_Budget_noBBP_extended_std, \
           Theoretical_Budget_Koestner,Theoretical_Budget_Koestner_std, \
           Theoretical_Budget_Koestner_extended,Theoretical_Budget_Koestner_extended_std, \
           POC_resp_mgC_m3_d_list,POC_resp_mgC_m3_d_std_list,bbp_PARR_dens0_densf,bbp_PARR_dens0_densf_std,\
           bbp_PARR_Koestner_dens0_densf,bbp_PARR_Koestner_dens0_densf_std,O2_resp_mgC_m3_d,O2_resp_mgC_m3_d_ci,list_Respi_types,n_profiles, \
           Delta_Integrated_POC, Delta_Integrated_POC_std, Delta_Integrated_POC_noBBP, Delta_Integrated_POC_noBBP_std, Delta_Integrated_POC_Koestner, Delta_Integrated_POC_Koestner_std, \
           Delta_flux, Delta_flux_std,Flux_dens0_mgC_m3_d,Flux_densf_mgC_m3_d, \
           depth_isopycnal,depth_isopycnal_down,depth_isopycnal_up,layer_thickness,MiP_POC_dens0_densf,MiP_POC_extended_dens0_densf,MaP_POC_dens0_densf,bbp_POC_Koestner_dens0_densf
# endregion
#######################################################################
# Parameters for the carbon budget calculation
#######################################################################
day0=datetime.datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime.datetime(2021,7,31)        # starting date for the carbon budget calculation
ndays=(dayf-day0).days          # number of days
dens00=1026.3                   # starting isopycnal
dens_thickness=0.05             # thickness of the layer considered (in kg/m3)
delta_dens=0.025                 # every time I do a loop, how much I do increase depth0
densff=1027.5                   # final isopycnal investigated

dens0_list=np.r_[dens00:densff-dens_thickness+0.01:delta_dens]

#######################################################################
# I loop on the different depths
#######################################################################
Theoretical_Budget_list = np.array([])
Theoretical_Budget_extended_list = np.array([])
Theoretical_Budget_std_list = np.array([])
Theoretical_Budget_extended_std_list = np.array([])
Theoretical_Budget_noBBP_list = np.array([])
Theoretical_Budget_noBBP_extended_list = np.array([])
Theoretical_Budget_noBBP_std_list = np.array([])
Theoretical_Budget_noBBP_extended_std_list = np.array([])
Theoretical_Budget_Koestner_list = np.array([])
Theoretical_Budget_Koestner_extended_list = np.array([])
Theoretical_Budget_Koestner_std_list = np.array([])
Theoretical_Budget_Koestner_extended_std_list = np.array([])
POC_resp_mgC_m3_d_list = np.array([])
POC_resp_mgC_m3_d_std_list = np.array([])
bbpPARR_mgC_m3_d_list = np.array([])
bbpPARR_mgC_m3_d_std_list = np.array([])
bbpPARR_Koestner_mgC_m3_d_list = np.array([])
bbpPARR_Koestner_mgC_m3_d_std_list = np.array([])
O2_resp_mgC_m3_d_list = np.array([])
O2_resp_mgC_m3_d_ci_list = np.array([])
depth_isopycnal_list = np.array([])
depth_isopycnal_down_list = np.array([])
depth_isopycnal_up_list = np.array([])
layer_thickness_list = np.array([])
dens0=dens0_list[0]
for dens0 in dens0_list:
    densf = dens0 + dens_thickness
    (Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_extended,Theoretical_Budget_extended_std,
       Theoretical_Budget_noBBP,Theoretical_Budget_noBBP_std,Theoretical_Budget_noBBP_extended,Theoretical_Budget_noBBP_extended_std,
       Theoretical_Budget_Koestner,Theoretical_Budget_Koestner_std,Theoretical_Budget_Koestner_extended,Theoretical_Budget_Koestner_extended_std,
       POC_resp_mgC_m3_d,POC_resp_mgC_m3_d_std,bbpPARR_mgC_m3_d,bbpPARR_mgC_m3_d_std,bbpPARR_Koestner_mgC_m3_d,bbpPARR_Koestner_mgC_m3_d_std,O2_resp_mgC_m3_d,O2_resp_mgC_m3_d_ci,RespirationTypes,n_profiles,
       _,_,_,_,_,_,_,_,_,_,depth_isopycnal,depth_isopycnal_down,depth_isopycnal_up,layer_thickness,_,_,_,_) = carbon_budget_calculation(dens0, densf, day0, dayf)

    Theoretical_Budget_list=np.append(Theoretical_Budget_list,Theoretical_Budget)
    Theoretical_Budget_extended_list=np.append(Theoretical_Budget_extended_list,Theoretical_Budget_extended)
    Theoretical_Budget_std_list=np.append(Theoretical_Budget_std_list,Theoretical_Budget_std)
    Theoretical_Budget_extended_std_list=np.append(Theoretical_Budget_extended_std_list,Theoretical_Budget_extended_std)
    Theoretical_Budget_noBBP_list=np.append(Theoretical_Budget_noBBP_list,Theoretical_Budget_noBBP)
    Theoretical_Budget_noBBP_extended_list=np.append(Theoretical_Budget_noBBP_extended_list,Theoretical_Budget_noBBP_extended)
    Theoretical_Budget_noBBP_std_list=np.append(Theoretical_Budget_noBBP_std_list,Theoretical_Budget_noBBP_std)
    Theoretical_Budget_noBBP_extended_std_list=np.append(Theoretical_Budget_noBBP_extended_std_list,Theoretical_Budget_noBBP_extended_std)
    Theoretical_Budget_Koestner_list=np.append(Theoretical_Budget_Koestner_list,Theoretical_Budget_Koestner)
    Theoretical_Budget_Koestner_extended_list=np.append(Theoretical_Budget_Koestner_extended_list,Theoretical_Budget_Koestner_extended)
    Theoretical_Budget_Koestner_std_list=np.append(Theoretical_Budget_Koestner_std_list,Theoretical_Budget_Koestner_std)
    Theoretical_Budget_Koestner_extended_std_list=np.append(Theoretical_Budget_Koestner_extended_std_list,Theoretical_Budget_Koestner_extended_std)
    POC_resp_mgC_m3_d_list=np.append(POC_resp_mgC_m3_d_list,POC_resp_mgC_m3_d,axis=0)
    POC_resp_mgC_m3_d_std_list=np.append(POC_resp_mgC_m3_d_std_list,POC_resp_mgC_m3_d_std,axis=0)
    bbpPARR_mgC_m3_d_list=np.append(bbpPARR_mgC_m3_d_list,bbpPARR_mgC_m3_d)
    bbpPARR_mgC_m3_d_std_list=np.append(bbpPARR_mgC_m3_d_std_list,bbpPARR_mgC_m3_d_std)
    bbpPARR_Koestner_mgC_m3_d_list=np.append(bbpPARR_Koestner_mgC_m3_d_list,bbpPARR_Koestner_mgC_m3_d)
    bbpPARR_Koestner_mgC_m3_d_std_list=np.append(bbpPARR_Koestner_mgC_m3_d_std_list,bbpPARR_Koestner_mgC_m3_d_std)
    O2_resp_mgC_m3_d_list=np.append(O2_resp_mgC_m3_d_list,O2_resp_mgC_m3_d)
    O2_resp_mgC_m3_d_ci_list=np.append(O2_resp_mgC_m3_d_ci_list,O2_resp_mgC_m3_d_ci.reshape((2,)),axis=0)
    depth_isopycnal_list=np.append(depth_isopycnal_list,depth_isopycnal)
    depth_isopycnal_down_list=np.append(depth_isopycnal_down_list,depth_isopycnal_down)
    depth_isopycnal_up_list=np.append(depth_isopycnal_up_list,depth_isopycnal_up)
    layer_thickness_list=np.append(layer_thickness_list,layer_thickness)

O2_resp_mgC_m3_d_ci_list=O2_resp_mgC_m3_d_ci_list.reshape(dens0_list.size,2)
POC_resp_mgC_m3_d_list=POC_resp_mgC_m3_d_list.reshape(dens0_list.size,len(RespirationTypes))
POC_resp_mgC_m3_d_std_list=POC_resp_mgC_m3_d_std_list.reshape(dens0_list.size,len(RespirationTypes))

dens_eddy_core_up = 1026.82
dens_eddy_core_down = 1027.2397618090454 #calculated at step 4 of Fig. 3a

#######################################################################
#region Plots
########################################################################################################################
######### Fig. 04A: with BBP from Koestner
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fs=9
width, height = 0.74, 0.4

fig, ax = plt.subplots(2, 1, figsize=(3.5, 3.5), sharex=True)
plt.subplots_adjust(hspace=0.05)
ax[0].set_position([0.18, 0.12+height+0.05, width, height])
ax[1].set_position([0.18 ,0.12, width, height])
ax[0].plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
ax[0].scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
ax[0].fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(2,3):
    ax[0].plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

ax[0].fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR')

ax[0].plot(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red')
ax[0].scatter(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red', s=5)
ax[0].fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_list - Theoretical_Budget_Koestner_std_list, Theoretical_Budget_Koestner_list + Theoretical_Budget_Koestner_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')
ax[0].hlines(200, xmin=ax[0].get_xlim()[0], xmax=ax[0].get_xlim()[1], color='darkgoldenrod')
ax[0].hlines(600, xmin=ax[0].get_xlim()[0], xmax=ax[0].get_xlim()[1], color='darkgoldenrod')
ax[0].hlines(depth_isopycnal_list[1], xmin=ax[0].get_xlim()[0], xmax=ax[0].get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
ax[0].set_xlim(-0.05,1)
ax[0].set_ylim(set_ylim_lower, set_ylim_upper)
ax[0].legend(fontsize=7)
ax[0].invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
    yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list)-1000 ))
ax[0].set_yticks(yticks)
ax[0].set_yticklabels(yticklabels,fontsize=6)
ax[0].text(0, .075, 'b', fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
# ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
ax[1].plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
ax[1].scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
ax[1].fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(9,10):
    ax[1].plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

ax[1].fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')

ax[1].plot(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red')
ax[1].scatter(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red', s=5)
ax[1].fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_extended_list - Theoretical_Budget_Koestner_extended_std_list, Theoretical_Budget_Koestner_extended_list + Theoretical_Budget_Koestner_extended_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')

ax[1].hlines(200, xmin=ax[1].get_xlim()[0], xmax=ax[1].get_xlim()[1], color='darkgoldenrod')
ax[1].hlines(600, xmin=ax[1].get_xlim()[0], xmax=ax[1].get_xlim()[1], color='darkgoldenrod')
ax[1].hlines(depth_isopycnal_list[1], xmin=ax[1].get_xlim()[0], xmax=ax[1].get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
ax[1].set_xlim(-0.05,1)
ax[1].set_ylim(set_ylim_lower, set_ylim_upper)
ax[1].set_xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
ax[1].invert_yaxis()
#I set yticks
ax[1].set_yticks(yticks)
ax[1].set_yticklabels(yticklabels,fontsize=6)
# ax.text(0.02, 1.075, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
ax[1].grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig04B_Koestner_v07.pdf' ,dpi=200)
plt.close()

########################################################################################################################
######### Supplementary Fig. 04A: with BBP from Koestner with Kalvealge Iversen and krem=0.1 d^-1
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(2,3):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')
plt.plot(POC_resp_mgC_m3_d_list[:, 0] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 5] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_list - Theoretical_Budget_Koestner_std_list, Theoretical_Budget_Koestner_list + Theoretical_Budget_Koestner_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
# plt.ylabel('Dens (kg/m$^3$)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
    yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig04A_Koestner_v07.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Supplementary Fig. 4B: Extended size spectrum with BBP from Koestner with Kalvealge Iversen and krem=0.1 d^-1
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fig = plt.figure(2, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(9,10):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m3_d_list[:, 7] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 12] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_extended_list - Theoretical_Budget_Koestner_extended_std_list, Theoretical_Budget_Koestner_extended_list + Theoretical_Budget_Koestner_extended_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')

plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
        yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(0.02, 1.075, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Fig04B_Koestner_v07.pdf' ,dpi=200)
plt.close()

########################################################################################################################
######### Supplementary Fig. 04A: with BBP from Cetinic
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 0], O2_resp_mgC_m3_d_ci_list[:, 1], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(2,3):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp]+bbpPARR_mgC_m3_d_list, depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_mgC_m3_d_list-POC_resp_mgC_m3_d_std_list[:,iResp]-bbpPARR_mgC_m3_d_std_list,
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_mgC_m3_d_list+POC_resp_mgC_m3_d_std_list[:,iResp]+bbpPARR_mgC_m3_d_std_list, facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')
plt.plot(POC_resp_mgC_m3_d_list[:, 0] + bbpPARR_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 5] + bbpPARR_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_list - Theoretical_Budget_std_list, Theoretical_Budget_list + Theoretical_Budget_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
# plt.ylabel('Dens (kg/m$^3$)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
    yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/Fig04A_Cetinic_v07.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Supplementary Fig. 04A: without BBP
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 0], O2_resp_mgC_m3_d_ci_list[:, 1], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(2,3):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp], depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp]-POC_resp_mgC_m3_d_std_list[:,iResp],
                  POC_resp_mgC_m3_d_list[:,iResp]+POC_resp_mgC_m3_d_std_list[:,iResp], facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')
plt.plot(POC_resp_mgC_m3_d_list[:, 0], depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 5], depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_noBBP_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_noBBP_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_noBBP_list - Theoretical_Budget_noBBP_std_list, Theoretical_Budget_noBBP_list + Theoretical_Budget_noBBP_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
# plt.ylabel('Dens (kg/m$^3$)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
    yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/Fig04A_noBBP_v07.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Supplementary Fig. Extended size spectrum: with BBP from Cetinic
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fig = plt.figure(2, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(9,10):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_mgC_m3_d_list , depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m3_d_list[:, 7] + bbpPARR_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 12] + bbpPARR_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_extended_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_extended_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_extended_list - Theoretical_Budget_extended_std_list, Theoretical_Budget_extended_list + Theoretical_Budget_extended_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')

plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
        yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/Suppl_CarbBudgetExtended_Cetinic_v07.pdf' ,dpi=200)
plt.close()

########################################################################################################################
######### Supplementary Fig. Extended size spectrum: without bbpPOC
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fig = plt.figure(2, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(9,10):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp], depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp]-POC_resp_mgC_m3_d_std_list[:,iResp],
                  POC_resp_mgC_m3_d_list[:,iResp]+POC_resp_mgC_m3_d_std_list[:,iResp], facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m3_d_list[:, 7], depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 12], depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_noBBP_extended_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_noBBP_extended_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_noBBP_extended_list - Theoretical_Budget_noBBP_extended_std_list, Theoretical_Budget_noBBP_extended_list + Theoretical_Budget_noBBP_extended_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')

plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,1)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
        yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v07/Supplementary/Suppl_CarbBudgetExtended_noBBP_v07.pdf' ,dpi=200)
plt.close()

#endregion
#######################################################################

#######################################################################
# I save some values for the latex document
#######################################################################
# region I save some values for the latex document
(Theoretical_Budget, Theoretical_Budget_std, Theoretical_Budget_extended, Theoretical_Budget_extended_std,Theoretical_Budget_noBBP,Theoretical_Budget_noBBP_std,
 Theoretical_Budget_noBBP_extended,Theoretical_Budget_noBBP_extended_std,Theoretical_Budget_Koestner, Theoretical_Budget_Koestner_std,
 Theoretical_Budget_Koestner_extended,Theoretical_Budget_Koestner_extended_std,POC_resp_mgC_m3_d, POC_resp_mgC_m3_d_std,
 bbpPARR_mgC_m3_d, bbpPARR_mgC_m3_d_std, bbpPARR_Koestner_mgC_m3_d,bbpPARR_Koestner_mgC_m3_d_std, O2_resp_mgC_m3_d, O2_resp_mgC_m3_d_ci, RespirationTypes,
 n_profiles, Delta_Integrated_POC, Delta_Integrated_POC_std,  Delta_Integrated_POC_noBBP, Delta_Integrated_POC_noBBP_std, Delta_Integrated_POC_Koestner,
 Delta_Integrated_POC_Koestner_std,Delta_flux, Delta_flux_std,Flux_dens0_mgC_m3_d,Flux_densf_mgC_m3_d, depth_isopycnal,
 depth_isopycnal_down, depth_isopycnal_up, layer_thickness,MiP_POC_dens0_densf,
 MiP_POC_extended_dens0_densf,MaP_POC_dens0_densf,bbp_POC_Koestner_dens0_densf) = carbon_budget_calculation(dens_eddy_core_up, dens_eddy_core_down, day0, dayf)

from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
#These values are extracted from the bulk POC and PARR calculated in the eddy core considering only one unique layer
argument = 'Flux_0413to0731_eddy_core_up'
arg_value=Flux_dens0_mgC_m3_d
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'Flux_0413to0731_eddy_core_down'
arg_value=Flux_densf_mgC_m3_d
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'Flux_0413to0731eddy_core_difference'
arg_value=abs((Flux_dens0_mgC_m3_d-Flux_densf_mgC_m3_d))
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'Flux_0413to0731eddy_core_difference_std'
arg_value=abs(Delta_flux_std)
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'Integrated_POC_0413to0731difference'
arg_value=abs(Delta_Integrated_POC_Koestner)
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'Integrated_POC_0413to0731difference_std'
arg_value=abs(Delta_Integrated_POC_Koestner_std)
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'bulkPOC_0413to0731_eddycore'
arg_value=Theoretical_Budget_Koestner
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'bulkPOC_0413to0731_eddycore_std'
arg_value=Theoretical_Budget_Koestner_std
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'bulkPOC_0413to0731_eddycore_extended'
arg_value=Theoretical_Budget_Koestner_extended
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'bulkPOC_0413to0731_eddycore_extended_std'
arg_value=Theoretical_Budget_Koestner_extended_std
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'PARR_0413to0731_eddycore'
arg_value=POC_resp_mgC_m3_d[2] + bbpPARR_Koestner_mgC_m3_d
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'PARR_0413to0731_eddycore_std'
arg_value=np.sqrt(POC_resp_mgC_m3_d_std[2]**2+bbpPARR_Koestner_mgC_m3_d_std**2)
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'PARR_0413to0731_eddycore_extended'
arg_value=POC_resp_mgC_m3_d[9] + bbpPARR_Koestner_mgC_m3_d
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'PARR_0413to0731_eddycore_extended_std'
arg_value=np.sqrt(POC_resp_mgC_m3_d_std[9]**2+bbpPARR_Koestner_mgC_m3_d_std**2)
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'PARR_0413to0731_eddycore_bbp_contribution'
arg_value=bbpPARR_Koestner_mgC_m3_d / (POC_resp_mgC_m3_d[2] + bbpPARR_Koestner_mgC_m3_d) *100
write_latex_data(filename,argument,'%d' % arg_value)
mipmap=np.nanmean(MiP_POC_dens0_densf)+np.nanmean(MaP_POC_dens0_densf)
mip_extended_only=np.nanmean(MiP_POC_extended_dens0_densf)-np.nanmean(MiP_POC_dens0_densf)
bbppoc=np.nanmean(bbp_POC_Koestner_dens0_densf)
argument = 'MiPextendedOnlyToMipMap_ratio_EC'
arg_value=np.round(mip_extended_only/mipmap)
write_latex_data(filename,argument,'%d' % arg_value)

# I calculate the oxygen cons rate as the mean of the oxygen cons rate values obtained in the eddy core
argument = 'OxyCR_0413to0731_eddycore'
# arg_value=O2_resp_mgC_m3_d
sel_eddycore = (depth_isopycnal_list>200)&(depth_isopycnal_list<600)
arg_value=np.mean(O2_resp_mgC_m3_d_list[sel_eddycore])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'OxyCR_0413to0731_eddycore_std'
# arg_value=abs(np.diff(O2_resp_mgC_m3_d_ci)[0][0]/2)
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list[sel_eddycore, 1] - O2_resp_mgC_m3_d_ci_list[sel_eddycore, 0])*0.5
write_latex_data(filename,argument,'%0.3f' % arg_value)

#These values are extracted from the bulk POC and PARR calculated for different layers within the eddy core
sel = (dens0_list+dens_thickness*0.5)>=dens_eddy_core_up
Theoretical_Budget_list_sel = Theoretical_Budget_list[sel]
Theoretical_Budget_std_list_sel = Theoretical_Budget_std_list[sel]
Theoretical_Budget_extended_list_sel = Theoretical_Budget_extended_list[sel]
Theoretical_Budget_extended_std_list_sel = Theoretical_Budget_extended_std_list[sel]
POC_resp_mgC_m3_d_list_sel = POC_resp_mgC_m3_d_list[sel,:]
POC_resp_mgC_m3_d_std_list_sel = POC_resp_mgC_m3_d_std_list[sel,:]
O2_resp_mgC_m3_d_list_sel = O2_resp_mgC_m3_d_list[sel]
O2_resp_mgC_m3_d_ci_list_sel = O2_resp_mgC_m3_d_ci_list[sel]
dens_list_sel = dens0_list[sel]+dens_thickness*0.5
depth_isopycnal_list_sel = depth_isopycnal_list[sel]
depth_isopycnal_down_list_sel = depth_isopycnal_down_list[sel]
depth_isopycnal_up_list_sel = depth_isopycnal_up_list[sel]
argument = 'bulkPOC_0413to0731_max_layer'
arg_value=Theoretical_Budget_list_sel.max()
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'bulkPOC_0413to0731_max_layer_std'
arg_value=Theoretical_Budget_std_list_sel[7]
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'bulkPOC_0413to0731_min_layer'
arg_value=Theoretical_Budget_list_sel.min()
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'bulkPOC_0413to0731_min_layer_std'
arg_value=Theoretical_Budget_std_list_sel[12]
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'PARR_0413to0731_max_layer'
arg_value=POC_resp_mgC_m3_d_list_sel[:,2].max()
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'PARR_0413to0731_max_layer_std'
arg_value=POC_resp_mgC_m3_d_std_list_sel[14,2]
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'PARR_0413to0731_min_layer'
arg_value=POC_resp_mgC_m3_d_list_sel[:,2].min()
write_latex_data(filename,argument,'%0.4f' % arg_value)
argument = 'PARR_0413to0731_min_layer_std'
arg_value=POC_resp_mgC_m3_d_std_list_sel[-1,2]
write_latex_data(filename,argument,'%0.4f' % arg_value)
print(depth_isopycnal_list_sel[0],depth_isopycnal_list_sel[8],dens_list_sel[0],dens_list_sel[8])
argument = 'OxyCR_0413to0731_175to295m_up'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[0:9,0])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'OxyCR_0413to0731_175to295m_down'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[0:9,1])
write_latex_data(filename,argument,'%0.3f' % arg_value)
print(depth_isopycnal_down_list_sel[10],depth_isopycnal_up_list_sel[10],dens_list_sel[10],dens_list_sel[10])
argument = 'OxyCR_0413to0731_320to365m_up'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[10,0])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'OxyCR_0413to0731_320to365m_down'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[10,1])
write_latex_data(filename,argument,'%0.3f' % arg_value)
print(depth_isopycnal_list_sel[11],depth_isopycnal_list_sel[14],dens_list_sel[11],dens_list_sel[14])
argument = 'OxyCR_0413to0731_373to484m_up'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[11:15,0])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'OxyCR_0413to0731_373to484m_down'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[11:15,1])
write_latex_data(filename,argument,'%0.3f' % arg_value)
print(depth_isopycnal_list_sel[15],depth_isopycnal_list_sel[17],dens_list_sel[15],dens_list_sel[17])
argument = 'OxyCR_0413to0731_520to581m_up'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[15:18,0])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'OxyCR_0413to0731_520to581m_down'
arg_value=np.mean(O2_resp_mgC_m3_d_ci_list_sel[15:18,1])
write_latex_data(filename,argument,'%0.3f' % arg_value)
print(depth_isopycnal_list_sel[11],depth_isopycnal_list_sel[12],dens_list_sel[11],dens_list_sel[12])
argument = 'bulkPOC_0413to0731_373to410m_extended'
arg_value=np.mean(Theoretical_Budget_extended_list_sel[11:13])
write_latex_data(filename,argument,'%0.3f' % arg_value)
argument = 'bulkPOC_0413to0731_373to410m_extended_std'
arg_value=np.mean(Theoretical_Budget_extended_std_list_sel[11:13])
write_latex_data(filename,argument,'%0.3f' % arg_value)

#These values are from literature
Oxy2C = 0.89                # to convert from mol of oxygen to mol of carbon
mol2gC = 12.0107            # to convert from mol of carbon to grams of carbon
argument = 'Kiko2020_zooRR_down'
arg_value=0.9*Oxy2C*mol2gC/1000
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Kiko2020_zooRR_up'
arg_value=2.8*Oxy2C*mol2gC/1000
write_latex_data(filename,argument,'%0.2f' % arg_value)
# endregion

# endregion

