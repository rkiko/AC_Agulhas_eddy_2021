import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
from datetime import date
import netCDF4 as nc
import seawater as sw
import gsw
from pathlib import Path
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import pandas as pd
import pickle
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

filename='6903095_Sprof.nc'

#######################################################################
# I load the data
#######################################################################

ds = nc.Dataset('%s/%s' % (storedir,filename))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])

Date_Num=np.array(ds.variables['JULD'])
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
date_reference_datenum=date.toordinal(date_reference)+366

Date_Vec=np.zeros([Date_Num.size,6])
for i in range(0,Date_Num.size):
    date_time_obj = date_reference + datetime.timedelta(days=Date_Num[i])
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

#######################################################################
# I calculate the Brunt Vaisala frequency, i.e. the square of the buoyancy frequency
#######################################################################
bvf=np.ones((temp.shape[0],temp.shape[1]))*99999
i=0
for i in range(0,bbp700.shape[0]):
    mask_bvf = np.logical_and(pres[i,:] != 99999, temp[i,:] != 99999, psal[i,:] != 99999)  # I exclude the points with value = 99999
    lat_tmp = np.tile(lat[i], [pres[i,:].shape[0], 1])
    lon_tmp = np.tile(lon[i], [pres[i,:].shape[0], 1])
    lat_tmp = np.squeeze(lat_tmp[mask_bvf])
    lon_tmp = np.squeeze(lon_tmp[mask_bvf])
    pres_tmp = pres[i,mask_bvf]
    psal_tmp = psal[i,mask_bvf]
    temp_tmp = temp[i,mask_bvf]
    abs_psal_tmp = gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp)  # I compute absolute salinity
    cons_tmp = gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)  # I compute conservative temperature
    bvf_tmp=gsw.Nsquared(abs_psal_tmp, cons_tmp, pres_tmp, lat=lat_tmp)[0]
    bvf[i,mask_bvf]= np.append(bvf_tmp,99999)

#######################################################################
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from paruvpy import bbp700toPOC
sPOC=bbp700.copy()*0+99999
i=0
for i in range(0,bbp700.shape[0]):
    bbp700tmp=bbp700[i,:]
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
    # I exclude nan values
    sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)
    bbp700tmp=bbp700tmp[sel]
    depth_tmp=depth_tmp[sel]
    temp_tmp=temp_tmp[sel]
    # I convert to small POC (sPOC) and I set to 0 values <0
    sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
    sPOC_tmp[sPOC_tmp<0]=0
    sPOC[i,sel]=sPOC_tmp

#######################################################################
# I load the eddy position calculated from satellite, which I need to calculate the distance with the BGC Argo
#######################################################################
filename_eddyData='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies.csv' % home
data_eddy=pd.read_csv(filename_eddyData, sep=',', header=0)
lonEddy=data_eddy['x_centroid']
latEddy=data_eddy['y_centroid']
Date_Num_Eddy=data_eddy['date_num']
radius_Vmax=data_eddy['Rmax']
del data_eddy
DateTime_Eddy=[]
i=Date_Num_Eddy[0]
for i in Date_Num_Eddy:
    DateTime_Eddy=np.append(DateTime_Eddy,datetime.datetime.fromordinal(int(i)-366))

#######################################################################
# I select the eddy position calculated from satellite in the days of the BGC Argo float profiles
#######################################################################

sel_Satel2Float=np.squeeze(np.zeros((1,Date_Num.size)))
i=0
for i in range(0,Date_Num.size):
    datetmp = np.squeeze(np.array([Date_Vec[i,0:3]]))
    datetmp = datetime.datetime(datetmp[0],datetmp[1],datetmp[2])
    idx=np.array(np.where(datetmp == DateTime_Eddy))
    if idx.size>1: print('error, i %d, ' % i, idx)
    if idx.size==0:
        print('Warning: missing eddy contour and center for %d-%d-%d' % (Date_Vec[i,0],Date_Vec[i,1],Date_Vec[i,2]))
        sel_Satel2Float[i] = i
        continue
    sel_Satel2Float[i] = idx[0]

sel_Satel2Float=sel_Satel2Float.astype(int)
sel_eddyCont_missing_days=sel_Satel2Float!=99999
lonEddy_4float=np.array(lonEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
latEddy_4float=np.array(latEddy[sel_Satel2Float[sel_eddyCont_missing_days]])
radius_Vmax_4float=np.array(radius_Vmax[sel_Satel2Float[sel_eddyCont_missing_days]])

#######################################################################
# I calculate the distance between the BGC argo and the eddy center
#######################################################################

dist=np.arccos(np.sin(lat*np.pi/180)*np.sin(latEddy_4float*np.pi/180)+np.cos(lat*np.pi/180)*np.cos(latEddy_4float*np.pi/180)*np.cos((lonEddy_4float-lon)*np.pi/180))*180/np.pi
dist_km=dist*111

#######################################################################
# I select the data only when the BGC Argo float was inside the eddy
#######################################################################
sel_insideEddy = dist_km <= radius_Vmax_4float

lon=lon[sel_insideEddy]
lat=lat[sel_insideEddy]
Date_Num=Date_Num[sel_insideEddy]
Date_Vec=Date_Vec[sel_insideEddy,:]
pres=pres[sel_insideEddy]
depth=depth[sel_insideEddy,:]
dens=dens[sel_insideEddy,:]
temp=temp[sel_insideEddy,:]
psal=psal[sel_insideEddy,:]
chla=chla[sel_insideEddy,:]
doxy=doxy[sel_insideEddy,:]
sPOC=sPOC[sel_insideEddy,:]
bvf=bvf[sel_insideEddy,:]

#######################################################################
# I load the euphotic layer depth
#######################################################################
a_file = open("%s/an36/data_an36.pkl" % storedir, "rb")
data_an36 = pickle.load(a_file)
zeu_float=data_an36['zeu_float']
zeu_datenum=data_an36['zeu_datenum']
zeu_datenum = zeu_datenum[~np.isnan(zeu_float)]
zeu_float = zeu_float[~np.isnan(zeu_float)]

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from paruvpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,chla.shape[0]):
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
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
critical_depth_1=data_an45['critical_depth_1']
critical_depth_2=data_an45['critical_depth_2']
critical_depth_datenum=data_an45['critical_depth_datenum']
critical_depth_datenum_1 = critical_depth_datenum[~np.isnan(critical_depth_1)]
critical_depth_datenum_2 = critical_depth_datenum[~np.isnan(critical_depth_2)]
critical_depth_datenum = critical_depth_datenum[~np.isnan(critical_depth)]
critical_depth = critical_depth[~np.isnan(critical_depth)]
critical_depth_1 = critical_depth_1[~np.isnan(critical_depth_1)]
critical_depth_2 = critical_depth_2[~np.isnan(critical_depth_2)]

#######################################################################
# I define the parameters list and I start the loop on each of them
#######################################################################
parameter=temp
parameter_ylabel_list=['Temperature ($^{\circ}$C)','Pratical salinity (psu)','Chlorophyll-a (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)','$N^2$ (s$^{-2}$)']
parameter_panellabel_list=['b','a','d','c','f',' ']
parameter_shortname_list=['temp','psal','chla','doxy','bbpPOC','BrVais']
ipar=3
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0: parameter=temp.copy()
    elif ipar==1:   parameter=psal.copy()
    elif ipar == 2: parameter=chla.copy()
    elif ipar == 3: parameter=doxy.copy()
    elif ipar == 4: parameter=sPOC.copy()
    elif ipar == 5: parameter=bvf.copy()

    #I filter the profiles
    parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([]);pressure_parameter=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        sel=(parameter[i,:]!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
        if sum(sel) > 0:
            z=parameter[i,sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
            z = savgol_filter(z, 5, 1)
            parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
            depth_parameter = np.concatenate((depth_parameter, y1));dens_parameter = np.concatenate((dens_parameter, y2));pressure_parameter = np.concatenate((pressure_parameter, y3))

    parameter_filtered[parameter_filtered<0]=0
    # I define the x and y arrays for the contourf plot
    x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
    y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
    y2_parameter = np.linspace(dens_parameter.min(),dens_parameter.max(),50)
    y3_parameter = np.linspace(pressure_parameter.min(),pressure_parameter.max(),50)
    # I interpolate
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
    parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y2_parameter)
    parameter_interp_dens = griddata((Date_Num_parameter,dens_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")


    ########################################################
    ####### I plot: versus depth
    ########################################################
    if ipar==4:
        parameter_interp_depth[parameter_interp_depth > 40] = 40
    elif ipar==5:
        parameter_interp_depth[parameter_interp_depth > 3*10**-5] = 3*10**-5

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
    plt.plot(critical_depth_datenum,critical_depth,'w');plt.plot(critical_depth_datenum,critical_depth,'w.')
    # plt.plot(critical_depth_datenum_1,critical_depth_1,'w--')#;plt.plot(critical_depth_datenum_1,critical_depth_1,'w.')
    # plt.plot(critical_depth_datenum_2,critical_depth_2,'w--')#;plt.plot(critical_depth_datenum_2,critical_depth_2,'w.')
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
    # I save
    plt.savefig('../Plots/an46/TimeSeries_%s_vs_01depth_an46.pdf' % parameter_shortname_list[ipar],dpi=200)
    plt.close()

#To plot SSI and K490
# fig = plt.figure(1, figsize=(12, 8))
#  ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, 500),
#                    xlim=(Date_Num.min(), Date_Num.max()))
# # ax_1 = plot2 = plt.contourf(x_parameter, y1_parameter, parameter_interp_depth)
# plt.plot(ssi_datenum, ssi_per_hour_float, 'b');
# plt.plot(ssi_datenum, Kd_490_float*1000, 'r');
# nxticks=10
# xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
# xticklabels=[]
# for i in xticks:
#     date_time_obj = date_reference + datetime.timedelta(days=i)
#     xticklabels.append(date_time_obj.strftime('%d %B'))
# ax.set_xticks(xticks)
# ax.set_xticklabels(xticklabels)
# plt.xticks(rotation=90,fontsize=12)
