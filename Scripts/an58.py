import numpy as np
import os
import netCDF4 as nc
import pickle
import seawater as sw
import calendar
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import gsw
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

########################################################################################################################
# I load the BGC Argo temperature, chl etc. profiles
########################################################################################################################
filename='6903095_Sprof.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
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

#######################################################################
#I tranform the pressure to depth
#######################################################################
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth=np.ones(pres.shape)*99999
depth[mask_depth]=depth_tmp

#######################################################################
#I compute the conservative temperature
#######################################################################
mask_dens = np.logical_and(pres != 99999, temp != 99999, psal != 99999)  # I exclude the points with value = 99999
lat_tmp = np.tile(lat, [pres.shape[1], 1]).T
lon_tmp = np.tile(lon, [pres.shape[1], 1]).T
lat_tmp = lat_tmp[mask_dens]
lon_tmp = lon_tmp[mask_dens]
pres_tmp = pres[mask_dens]
psal_tmp = psal[mask_dens]
temp_tmp = temp[mask_dens]
abs_psal_tmp = gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp)  # I compute absolute salinity
cons_tmp = gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)  # I compute conservative temperature
abs_psal=np.ones(temp.shape)*99999
abs_psal[mask_dens]=abs_psal_tmp
cons_temp=np.ones(temp.shape)*99999
cons_temp[mask_dens]=cons_tmp

########################################################################################################################
# I load the eddy radiuses and the eddy-float distance
########################################################################################################################
a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an57/data_an57.pkl" % (home), "rb")
data_an57 = pickle.load(a_file)
Date_Num_float = data_an57['Date_Num_float']
dist_float_eddy_km = data_an57['dist_float_eddy_km']
radius_Vmax_floatDays = data_an57['radius_Vmax_floatDays']
radius_Out_floatDays = data_an57['radius_Out_floatDays']
a_file.close()

########################################################################################################################
# I plot by doing a loop on the different profiles
########################################################################################################################
fs=10
width, height = 0.73, 0.75
set_ylim_lower, set_ylim_upper = 0, 600
temp0, temp1 = 4,15
doxy0, doxy1 = 180,245

fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.gca().invert_yaxis()
fig = plt.figure(2, figsize=(3.5, 3.5))
ax2 = fig.add_axes([0.18, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.gca().invert_yaxis()
fig = plt.figure(3, figsize=(3.5, 3.5))
ax2 = fig.add_axes([0.18, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.gca().invert_yaxis()

i=0
for i in range (0,temp.shape[0]):
    dist_float_eddy_km_tmp = dist_float_eddy_km[i]
    radius_Vmax_floatDays_tmp = radius_Vmax_floatDays[i]
    if dist_float_eddy_km_tmp>radius_Vmax_floatDays_tmp:    continue

    temp_tmp=temp[i,:]
    cons_temp_tmp=cons_temp[i,:]
    doxy_tmp=doxy[i,:]
    depth_tmp=depth[i,:]
    sel = (depth_tmp != 99999) & (temp_tmp != 99999) & (cons_temp_tmp != 99999) & (doxy_tmp != 99999)
    depth_tmp = depth_tmp[sel]
    temp_tmp = temp_tmp[sel]
    cons_temp_tmp = cons_temp_tmp[sel]
    doxy_tmp = doxy_tmp[sel]

    dist_float_eddy_km_tmp=np.tile(dist_float_eddy_km_tmp,temp_tmp.size)
    temp_tmp[temp_tmp<temp0]=temp0
    temp_tmp[temp_tmp>temp1]=temp1
    cons_temp_tmp[cons_temp_tmp<temp0]=temp0
    cons_temp_tmp[cons_temp_tmp>temp1]=temp1
    doxy_tmp[doxy_tmp<doxy0]=doxy0
    doxy_tmp[doxy_tmp>doxy1]=doxy1
    ####################################################################################################################
    # Plot part
    ####################################################################################################################
    plt.figure(1)
    plot1 = plt.scatter(dist_float_eddy_km_tmp, depth_tmp, c=temp_tmp,s=5)
    plt.figure(2)
    plot2 = plt.scatter(dist_float_eddy_km_tmp, depth_tmp, c=doxy_tmp,s=5)
    plt.figure(3)
    plot2 = plt.scatter(dist_float_eddy_km_tmp, depth_tmp, c=cons_temp_tmp,s=5)


# I add last details
plt.figure(1)
plt.xlabel('Distance from center (km)', fontsize=fs)
plt.ylabel('Depth (m)', fontsize=fs)
cbar = plt.colorbar(plot1)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Temp. (°C)', fontsize=fs)
plt.savefig('../Plots/an58/01Temp_profiles_vs_distance_from_center_an58.pdf' ,dpi=200)
plt.close()
plt.figure(2)
plt.xlabel('Distance from center (km)', fontsize=fs)
plt.ylabel('Depth (m)', fontsize=fs)
cbar = plt.colorbar(plot2)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Doxy ($\mu$mol/kg)', fontsize=fs)
plt.savefig('../Plots/an58/02Doxy_profiles_vs_distance_from_center_an58.pdf' ,dpi=200)
plt.close()
plt.figure(3)
plt.xlabel('Distance from center (km)', fontsize=fs)
plt.ylabel('Depth (m)', fontsize=fs)
cbar = plt.colorbar(plot1)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Cons. Temp. (°C)', fontsize=fs)
plt.savefig('../Plots/an58/03ConsTemp_profiles_vs_distance_from_center_an58.pdf' ,dpi=200)
plt.close()






