import os
import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
import random
import gsw
from pathlib import Path
home = str(Path.home())
#globals().clear()
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

filename='6903095_Sprof.nc'

#######################################################################
# I load the data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
Date_Num_float = Date_Num + matlab_datenum(1950,1,1)
pres = np.array(ds.variables['PRES'])
psal = np.array(ds.variables['PSAL'])
temp = np.array(ds.variables['TEMP'])

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

#######################################################################
# I plot TEMP vs. PSAL for each profile
#######################################################################

list_distances=[radius_Vmax_floatDays.max(),40,30,20]
width, height = 0.65, 0.65

for distance_limit in list_distances:
    fig = plt.figure(1, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))
    fig = plt.figure(2, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))

    i=0
    for i in range(0, pres.shape[0]):
        dist_float_eddy_km_tmp = dist_float_eddy_km[i]
        if dist_float_eddy_km_tmp > distance_limit:    continue

        x = psal[i,:]
        y = temp[i,:]
        sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
        x=x[sel];y=y[sel]
        Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
        plt.figure(1)
        plot1=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])
        x = abs_psal[i,:]
        y = cons_temp[i,:]
        sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
        x=x[sel];y=y[sel]
        Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
        plt.figure(2)
        plot2=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])

    plt.figure(1)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.xlabel('Salinity (PSU)', fontsize=10)
    plt.ylabel('Temperature (°C)', fontsize=10)
    plt.title('Profiles less than\n%d km from eddy center' % distance_limit, fontsize=10)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an61/01TempVsPsal_%03dkm_an61.pdf' % distance_limit, dpi=200)
    plt.close()

    plt.figure(2)
    cbar = plt.colorbar(plot2)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.xlabel('Abs. Salinity (PSU)', fontsize=10)
    plt.ylabel('Cons. Temperature (°C)', fontsize=10)
    plt.title('Profiles less than\n%d km from eddy center' % distance_limit, fontsize=10)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an61/02ConsTempVsAbsPsal_%03dkm_an61.pdf' % distance_limit, dpi=200)
    plt.close()

    #Profiles in random order
    a=np.array(range(0, pres.shape[0]))
    random.shuffle(a)

    fig = plt.figure(1, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))
    fig = plt.figure(2, figsize=(3, 3))
    ax2 = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))
    i=0
    for i in a:
        dist_float_eddy_km_tmp = dist_float_eddy_km[i]
        radius_Vmax_floatDays_tmp = radius_Vmax_floatDays[i]
        if dist_float_eddy_km_tmp > radius_Vmax_floatDays_tmp:    continue

        x = psal[i,:]
        y = temp[i,:]
        sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
        x=x[sel];y=y[sel]
        Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
        plt.figure(1)
        plot1=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])
        x = abs_psal[i,:]
        y = cons_temp[i,:]
        sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
        x=x[sel];y=y[sel]
        Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
        plt.figure(2)
        plot2=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])

    plt.figure(1)
    cbar = plt.colorbar(plot1)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.xlabel('Salinity (PSU)', fontsize=10)
    plt.ylabel('Temperature (°C)', fontsize=10)
    plt.title('Profiles less than\n%d km from eddy center' % distance_limit, fontsize=10)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an61/01TempVsPsal_random_%03dkm_an61.pdf' % distance_limit, dpi=200)
    plt.close()

    plt.figure(2)
    cbar = plt.colorbar(plot2)
    nxticks = 10
    xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
    cbar.set_ticks(xticks)
    cbar.set_ticklabels(xticklabels)
    plt.xlabel('Abs. Salinity (PSU)', fontsize=10)
    plt.ylabel('Cons. Temperature (°C)', fontsize=10)
    plt.title('Profiles less than\n%d km from eddy center' % distance_limit, fontsize=10)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an61/02ConsTempVsAbsPsal_random_%03dkm_an61.pdf' % distance_limit, dpi=200)
    plt.close()


