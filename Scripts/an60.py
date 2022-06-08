import numpy as np
import os
import netCDF4 as nc
import pickle
import seawater as sw
import gsw
import calendar
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
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
# I plot by doing two loops: on on the different layers and one on the different profiles
########################################################################################################################

# General parameters
depth_tmp0_list =np.array([100,200,300,400,500,600])
delta_depth = 15 #how many meters I consider around each layer
# Parameters for the plot
fs=10
width, height = 0.73, 0.75
set_xlim_lower, set_xlim_upper = 0, 50
temp0, temp1 = 3,16
doxy0, doxy1 = 180,245

idepth=0
for idepth in range(0,depth_tmp0_list.size):
    depth_tmp0=depth_tmp0_list[idepth]

    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height], ylim=(temp0, temp1),xlim=(set_xlim_lower, set_xlim_upper))
    fig = plt.figure(2, figsize=(3.5, 3.5))
    ax2 = fig.add_axes([0.18, 0.15, width, height], ylim=(doxy0, doxy1),xlim=(set_xlim_lower, set_xlim_upper))
    fig = plt.figure(3, figsize=(3.5, 3.5))
    ax3 = fig.add_axes([0.18, 0.15, width, height], ylim=(temp0, temp1),xlim=(set_xlim_lower, set_xlim_upper))

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

        sel_depth = (depth_tmp >=(depth_tmp0 - delta_depth))&(depth_tmp <(depth_tmp0 + delta_depth))
        if sum(sel_depth)==0: continue
        temp_tmp = np.mean(temp_tmp[sel_depth])
        cons_temp_tmp = np.mean(cons_temp_tmp[sel_depth])
        doxy_tmp = np.mean(doxy_tmp[sel_depth])

        Date_Num_float_tmp=Date_Num_float[i]
        ####################################################################################################################
        # Plot part
        ####################################################################################################################
        plt.figure(1)
        plot1 = plt.scatter(dist_float_eddy_km_tmp, temp_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])
        plt.figure(2)
        plot2 = plt.scatter(dist_float_eddy_km_tmp, doxy_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])
        plt.figure(3)
        plot3 = plt.scatter(dist_float_eddy_km_tmp, cons_temp_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])


    # I add last details
    plt.figure(1)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Temp. (°C)', fontsize=fs)
    plt.title('%d m [%d—%d]' % (depth_tmp0,depth_tmp0-delta_depth,depth_tmp0+delta_depth), fontsize=fs)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('../Plots/an60/01Temp_profiles_vs_distance_from_center_%dm_an60.pdf' % (depth_tmp0) ,dpi=200)
    plt.close()

    plt.figure(2)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Doxy ($\mu$mol/kg)', fontsize=fs)
    plt.title('%d m [%d—%d]' % (depth_tmp0,depth_tmp0-delta_depth,depth_tmp0+delta_depth), fontsize=fs)
    cbar = plt.colorbar(plot2)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('../Plots/an60/02Doxy_profiles_vs_distance_from_center_%dm_an60.pdf' % (depth_tmp0) ,dpi=200)
    plt.close()

    plt.figure(3)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Cons. temp. (°C)', fontsize=fs)
    plt.title('%d m [%d—%d]' % (depth_tmp0,depth_tmp0-delta_depth,depth_tmp0+delta_depth), fontsize=fs)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('../Plots/an60/03ConsTemp_profiles_vs_distance_from_center_%dm_an60.pdf' % (depth_tmp0) ,dpi=200)
    plt.close()






