import os,sys
import numpy as np
import netCDF4 as nc
import pickle
import seawater as sw
import gsw
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import datetime
from pathlib import Path
home = str(Path.home())
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

filename='6903095_Sprof_all.nc'

#######################################################################
# I load the data
#######################################################################

ds = nc.Dataset('%s/%s' % (storedir,filename))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])

Date_Num=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)

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

#######################################################################
# I compute the potential density: for that, I need absolute salinity and
# conservative temperature, so I transform salinity and temperature first
#######################################################################
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
cons_temp=np.ones(temp.shape)*99999
cons_temp[mask_dens]=cons_tmp
dens=np.ones(temp.shape)*99999
dens[mask_dens]=dens_tmp+1000

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from oceanpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,temp.shape[0]):
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
    # I exclude nan values
    sel_non_nan=(depth_tmp!=99999)&(temp_tmp!=99999)
    temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
    mld_tmp,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
    mld=np.append(mld,mld_tmp)

#######################################################################
# For each profile, I calculate the isopycnal associated with the mixed layer depth
#######################################################################
isopycnal_mld=np.zeros((temp.shape[0],))
i=0
for i in range(0,temp.shape[0]):
    depth_tmp = depth[i,:]
    dens_tmp = dens[i, :]
    mld_tmp = mld[i]
    a=depth_tmp-mld_tmp
    idx=np.sum(a<0)-1 #index corresponding to the depth just slightly shallower than the MLD
    dist=(mld_tmp-depth_tmp[idx])/(depth_tmp[idx+1]-depth_tmp[idx])
    isopycnal_mld[i] = dens_tmp[idx] + (dens_tmp[idx+1] - dens_tmp[idx])*dist

#######################################################################
# For a given series of isopycnals, I calculate the depth at which they are found
#######################################################################
delta_isopycnal=0.05
delta_isopycnal2=0.05
list_isopycnals=np.r_[1026:1027+delta_isopycnal*0.5:delta_isopycnal]

isopycnal_vs_depth=np.zeros((list_isopycnals.size,cons_temp.shape[0]))
i=0
for i in range(0,list_isopycnals.size):
    isopycnal_tmp=list_isopycnals[i]
    isopycnal_tmp0=list_isopycnals[i]-delta_isopycnal2*0.5
    isopycnal_tmp1=list_isopycnals[i]+delta_isopycnal2*0.5
    # if i==0:    isopycnal_tmp0=isopycnal_tmp
    # if i==(list_isopycnals.size-1): isopycnal_tmp1=isopycnal_tmp
    j=0
    for j in range(0,cons_temp.shape[0]):
        depth_tmp = depth[j, :]
        dens_tmp = dens[j, :]
        sel_isopycnal = (dens_tmp>=isopycnal_tmp0)&(dens_tmp<isopycnal_tmp1)
        if sum(sel_isopycnal)==0:
            index = -1
            k=0
            for k in range(0,cons_temp.shape[1]-1):
                dens_0=dens_tmp[k]
                dens_1=dens_tmp[k+1]
                if (dens_0<isopycnal_tmp)&(dens_1>=isopycnal_tmp):
                    index=k;break

            if index==-1:
                if isopycnal_tmp<dens_tmp.min():
                    isopycnal_vs_depth[i, j] = 0
                else:
                    isopycnal_vs_depth[i, j] = depth_tmp.max()
            else:
                dist=(isopycnal_tmp-dens_0)/(dens_1-dens_0)
                isopycnal_vs_depth[i, j] = depth_tmp[k]+(depth_tmp[k+1]-depth_tmp[k])*dist

        else:
            isopycnal_vs_depth[i,j] = np.mean(depth_tmp [sel_isopycnal])

#######################################################################
# I save: I manually assign the isopycnal of the mixed layer depth based on the first of the two plots above
#######################################################################
dictionary_data = {"mld": mld,"isopycnal_mld": 1026.35,"isopycnal_mld_vs_time": isopycnal_mld,"Date_Num": Date_Num,"lon": lon,"lat": lat}
a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an68/data_MLD_an68.pkl" % (home), "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()

#######################################################################
# Plot part: I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)

filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

cons_temp=cons_temp[sel_insideEddy,:]
Date_Num=Date_Num[sel_insideEddy]
depth=depth[sel_insideEddy,:]
dens=dens[sel_insideEddy,:]
mld=mld[sel_insideEddy]
isopycnal_vs_depth=isopycnal_vs_depth[:,sel_insideEddy]
isopycnal_mld=isopycnal_mld[sel_insideEddy]

isopycnal_mld_mean=np.mean(isopycnal_mld)

#######################################################################
# Plot part: I interopolate the conservative temperature
#######################################################################
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")

parameter_ylabel_list=['Temperature ($^{\circ}$C)']
parameter_panellabel_list=['b']
parameter_shortname_list=['cons_temp']

ipar=0
parameter=cons_temp.copy()
#I filter the profiles
parameter_filtered=np.array([]);dens_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([])
i=0
for i in range(0,parameter.shape[0]):
    z = parameter[i, :]
    sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
    if ipar == 3: sel = (sel) & (z <= 100)
    if sum(sel) > 0:
        z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
        z = savgol_filter(z, 5, 1)
        parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
        depth_parameter = np.concatenate((depth_parameter, y1));dens_parameter = np.concatenate((dens_parameter, y2));
        y2 = savgol_filter(y2, 5, 1)
        dens_filtered = np.concatenate((dens_filtered, y2))

parameter_filtered[parameter_filtered<0]=0
# dens_filtered[dens_filtered<1026]=1026
dens_filtered[dens_filtered>1027.5]=1027.5
# I define the x and y arrays for the contourf plot
x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
# I interpolate
x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
dens_interp_depth = griddata((Date_Num_parameter,depth_parameter), dens_filtered, (x_parameter_g, y_parameter_g), method="nearest")
parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")


########################################################
####### I plot: versus depth
########################################################
import matplotlib.pyplot as plt

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
plt.plot(Date_Num,mld,'w',label='MLD')
i=0;plt.plot(Date_Num,isopycnal_vs_depth[i,:],label='%0.1f' % list_isopycnals[i])
i=6;plt.plot(Date_Num,isopycnal_vs_depth[i,:],label='%0.1f' % list_isopycnals[i])
i=7;plt.plot(Date_Num,isopycnal_vs_depth[i,:],label='%0.2f' % list_isopycnals[i])
# plot3 = ax.contour(x_parameter,y1_parameter, dens_interp_depth,levels=np.array([isopycnal_mld_mean,1026.3,1026.35]),colors='black',linestyles='solid',linewidths=1,zorder=10)#,cmap='Blues_r')
# plot3 = ax.contour(x_parameter,y1_parameter, dens_interp_depth,levels=np.r_[isopycnal_mld_mean:1026.4:(1026.4-isopycnal_mld_mean-0.000001)],colors='black',linestyles='solid',linewidths=1,zorder=10)#,cmap='Blues_r')
# plot4 = ax.contour(x_parameter,y1_parameter, dens_interp_depth,levels=np.r_[1026.4],colors='black',linestyles='dashed',linewidths=1,zorder=11)#,cmap='Blues_r')
# ax.clabel(plot3, inline=1, fontsize=10)
# ax.clabel(plot4, inline=1, fontsize=10)

plt.gca().invert_yaxis()
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
plt.legend(fontsize=15)
# I save
plt.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an68/Temp_density_an68.pdf' % home,dpi=200)
plt.close()


fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height])#, ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(Date_Num,isopycnal_mld,'.b-')
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
plt.ylabel('Density (kg/m3)')
plt.title('Isopycnal at the mixed layer depth VS time')
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# I save
plt.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an68/Isopycnal_of_MLD_vs_time_an68.pdf' % home,dpi=200)
plt.close()

