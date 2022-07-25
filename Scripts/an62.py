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
filename='6903095_Sprof_all.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
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
#I compute the conservative temperature, the absolute salinity, and the density
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
dens_tmp = gsw.density.sigma0(abs_psal_tmp, cons_tmp)
abs_psal=np.ones(temp.shape)*99999
abs_psal[mask_dens]=abs_psal_tmp
cons_temp=np.ones(temp.shape)*99999
cons_temp[mask_dens]=cons_tmp
dens = np.ones(temp.shape) * 99999
dens[mask_dens] = dens_tmp + 1000

########################################################################################################################
# I load the eddy radiuses and the eddy-float distance
########################################################################################################################
day_end_timeseries=np.array([2021,8,1])
image_path='../Plots/an62/20210413to%d%02d%02d' % (day_end_timeseries[0],day_end_timeseries[1],day_end_timeseries[2])
if not os.path.isdir(image_path):   os.system('mkdir %s' % image_path)
day_end_timeseries=matlab_datenum(day_end_timeseries)
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)
dist_float_eddy_km = data_dist_radius['Distance_Centroid']
sel_insideEddy = data_dist_radius['sel_insideEddy']
sel_insideEddy = (Date_Num<=day_end_timeseries)&(sel_insideEddy==1)
del data_dist_radius


########################################################################################################################
# I plot by doing two loops: one on the different isopycnals and one on the different profiles
########################################################################################################################

# General parameters
delta_rho = 0.05
reference_isopycnal_list = np.r_[1026.8:1027.50001:delta_rho]

depth_tmp0_list =np.array([100,200,300,400,500,600])
delta_depth = 15 #how many meters I consider around each layer
# Parameters for the plot
fs=10
width, height = 0.73, 0.75
set_xlim_lower, set_xlim_upper = 0, 50
temp0, temp1 = 3,8
doxy0, doxy1 = 180,245

i_isop=0
for i_isop in range(0,reference_isopycnal_list.size):
    reference_isopycnal = reference_isopycnal_list[i_isop]
    reference_isopycnal_down = reference_isopycnal - delta_rho / 2
    if reference_isopycnal_down < reference_isopycnal_list[0]:  reference_isopycnal_down = reference_isopycnal_list[0]
    reference_isopycnal_up = reference_isopycnal + delta_rho / 2
    if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]

    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height], ylim=(temp0, temp1),xlim=(set_xlim_lower, set_xlim_upper))
    fig = plt.figure(2, figsize=(3.5, 3.5))
    ax2 = fig.add_axes([0.18, 0.15, width, height], ylim=(doxy0, doxy1),xlim=(set_xlim_lower, set_xlim_upper))
    fig = plt.figure(3, figsize=(3.5, 3.5))
    ax3 = fig.add_axes([0.18, 0.15, width, height], ylim=(temp0, temp1),xlim=(set_xlim_lower, set_xlim_upper))

    max_temp=-99999;min_temp=99999
    max_cons_temp=-99999;min_cons_temp=99999
    max_doxy=-99999;min_doxy=99999
    depth_isopycnal_tmp = np.array([])
    i=0
    for i in range (0,Date_Num.size):
        dist_float_eddy_km_tmp = dist_float_eddy_km[i]
        if not sel_insideEddy[i]:    continue

        temp_tmp=temp[i,:]
        cons_temp_tmp=cons_temp[i,:]
        doxy_tmp=doxy[i,:]
        depth_tmp=depth[i,:]
        dens_tmp=dens[i,:]
        sel = (depth_tmp != 99999) & (dens_tmp != 99999) & (temp_tmp != 99999) & (cons_temp_tmp != 99999) & (doxy_tmp != 99999)
        depth_tmp = depth_tmp[sel]
        dens_tmp = dens_tmp[sel]
        temp_tmp = temp_tmp[sel]
        cons_temp_tmp = cons_temp_tmp[sel]
        doxy_tmp = doxy_tmp[sel]

        sel_layer = (dens_tmp >= reference_isopycnal_down) & (dens_tmp < reference_isopycnal_up)
        if sum(sel_layer)==0: continue
        temp_tmp = np.mean(temp_tmp[sel_layer])
        cons_temp_tmp = np.mean(cons_temp_tmp[sel_layer])
        doxy_tmp = np.mean(doxy_tmp[sel_layer])
        depth_tmp = np.mean(depth_tmp[sel_layer])
        depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, depth_tmp)

        max_temp=np.max( np.array( [max_temp,temp_tmp] ) );min_temp=np.min( np.array( [min_temp,temp_tmp] ) )
        max_cons_temp=np.max( np.array( [max_cons_temp,cons_temp_tmp] ) );min_cons_temp=np.min( np.array( [min_cons_temp,cons_temp_tmp] ) )
        max_doxy=np.max( np.array( [max_doxy,doxy_tmp] ) );min_doxy=np.min( np.array( [min_doxy,doxy_tmp] ) )

        Date_Num_float_tmp=Date_Num[i]
        ####################################################################################################################
        # Plot part
        ####################################################################################################################
        plt.figure(1)
        plot1 = plt.scatter(dist_float_eddy_km_tmp, temp_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num.min(),vmax=day_end_timeseries)
        plt.figure(2)
        plot2 = plt.scatter(dist_float_eddy_km_tmp, doxy_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num.min(),vmax=day_end_timeseries)
        plt.figure(3)
        plot3 = plt.scatter(dist_float_eddy_km_tmp, cons_temp_tmp, c=Date_Num_float_tmp,s=5,vmin=Date_Num.min(),vmax=day_end_timeseries)


    # I add last details
    plt.figure(1)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Temp. (°C)', fontsize=fs)
    plt.title('%d±%d m [%0.2f kg/m$^3$]' % (depth_isopycnal_tmp.mean(),np.round(depth_isopycnal_tmp.std()),reference_isopycnal), fontsize=fs)
    tmp=(max_temp-min_temp)/10
    plt.ylim(min_temp-tmp,max_temp+tmp)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num.min(), day_end_timeseries, nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('%s/01Temp_profiles_vs_distance_from_center_%0.2fkgm3_an62.pdf' % (image_path,reference_isopycnal) ,dpi=200)
    plt.close()

    plt.figure(2)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Doxy ($\mu$mol/kg)', fontsize=fs)
    plt.title('%d±%d m [%0.2f kg/m$^3$]' % (depth_isopycnal_tmp.mean(),np.round(depth_isopycnal_tmp.std()),reference_isopycnal), fontsize=fs)
    tmp=(max_doxy-min_doxy)/10
    plt.ylim(min_doxy-tmp,max_doxy+tmp)
    cbar = plt.colorbar(plot2)
    nxticks = 10
    xticks = np.linspace(Date_Num.min(), day_end_timeseries, nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('%s/02Doxy_profiles_vs_distance_from_center_%0.2fkgm3_an62.pdf' % (image_path,reference_isopycnal) ,dpi=200)
    plt.close()

    plt.figure(3)
    plt.xlabel('Distance from center (km)', fontsize=fs)
    plt.ylabel('Cons. temp. (°C)', fontsize=fs)
    plt.title('%d±%d m [%0.2f kg/m$^3$]' % (depth_isopycnal_tmp.mean(),np.round(depth_isopycnal_tmp.std()),reference_isopycnal), fontsize=fs)
    tmp=(max_cons_temp-min_cons_temp)/10
    plt.ylim(min_cons_temp-tmp,max_cons_temp+tmp)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num.min(), day_end_timeseries, nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2],matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.savefig('%s/03ConsTemp_profiles_vs_distance_from_center_%0.2fkgm3_an62.pdf' % (image_path,reference_isopycnal) ,dpi=200)
    plt.close()






