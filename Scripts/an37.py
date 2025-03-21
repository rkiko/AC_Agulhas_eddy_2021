import os
import ftplib
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
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from O2sat_func import O2sat_func
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

filename='6903095_Sprof_all.nc'

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
from oceanpy import bbp700toPOC
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
cons_temp=cons_temp[sel_insideEddy,:]
abs_psal=abs_psal[sel_insideEddy,:]
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
from oceanpy import mixed_layer_depth
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
# I define the parameters list and I start the loop on each of them
#######################################################################
parameter=temp
parameter_ylabel_list=['Temperature ($^{\circ}$C)','Temperature ($^{\circ}$C)','Pratical salinity (psu)','Absolute salinity (g/kg)','Chlorophyll-a (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)','$N^2$ (s$^{-2}$)']
parameter_panellabel_list=['b','b','a','a','c','d','f',' ']
parameter_shortname_list=['temp','cons_temp','psal','abs_psal','chla','doxy','bbpPOC','BrVais']
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)
ipar=7
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0: parameter=temp.copy()
    elif ipar==1:   parameter=cons_temp.copy()
    elif ipar==2:   parameter=psal.copy()
    elif ipar==3:   parameter=abs_psal.copy()
    elif ipar == 4: parameter=chla.copy()
    elif ipar == 5: parameter=doxy.copy()
    elif ipar == 6: parameter=sPOC.copy()
    elif ipar == 7: parameter=bvf.copy()

    #I filter the profiles
    parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        sel=(parameter[i,:]!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
        if sum(sel) > 0:
            z=parameter[i,sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
            z = savgol_filter(z, 5, 1)
            parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
            depth_parameter = np.concatenate((depth_parameter, y1));dens_parameter = np.concatenate((dens_parameter, y2))

    parameter_filtered[parameter_filtered<0]=0
    # I define the x and y arrays for the contourf plot
    x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
    y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
    y2_parameter = np.linspace(dens_parameter.min(),dens_parameter.max(),50)
    # I interpolate
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
    parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
    dens_interp_depth = griddata((Date_Num_parameter,depth_parameter), dens_parameter, (x_parameter_g, y_parameter_g), method="nearest")
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y2_parameter)
    parameter_interp_dens = griddata((Date_Num_parameter,dens_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")

    if ipar==0: temp_interp = parameter_interp_depth.copy()
    elif ipar==1:   cons_temp_interp_depth = parameter_interp_depth.copy();cons_temp_interp_dens = parameter_interp_dens.copy()
    elif ipar==2:   psal_interp = parameter_interp_depth.copy()
    elif ipar==3:   abs_psal_interp = parameter_interp_depth.copy()
    elif ipar == 4: chla_interp = parameter_interp_depth.copy()
    elif ipar == 5: doxy_interp = parameter_interp_depth.copy()
    elif ipar == 6: sPOC_interp = parameter_interp_depth.copy()
    elif ipar == 7: bvf_interp = parameter_interp_depth.copy()

    ########################################################
    ####### I plot: versus depth
    ########################################################
    if ipar==6:
        parameter_interp_depth[parameter_interp_depth > 40] = 40
    elif ipar==7:
        parameter_interp_depth[parameter_interp_depth > 3*10**-5] = 3*10**-5

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
    if (ipar==0)|(ipar==1):
        plt.plot(Date_Num,mld,'w')
    elif ipar==4:
        plt.plot(zeu_datenum,zeu_float,'w')
    elif ipar==7:
        plt.plot(Date_Num, mld, 'orangered')
        plot3 = ax.contour(x_parameter, y1_parameter, dens_interp_depth, levels=[1026.35],colors='orangered', linestyles='dashed', linewidths=1, zorder=10 )#,cmap='RdBu')
        fmt = {}
        strs = ['1026.35 kg/m$^3$']
        for l, s in zip(plot3.levels, strs):
            fmt[l] = s
        ax.clabel(plot3, plot3.levels[::], inline=True, fmt=fmt, fontsize=10)

    plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w')
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
    plt.savefig('../Plots/an37/TimeSeries_%s_vs_01depth_an37.pdf' % parameter_shortname_list[ipar],dpi=200)
    plt.close()


    #Here, for the temperature and oxygen, I analyse how much they changed in the top 200 and in the 200-600 m layers
    if (ipar == 1)|(ipar == 4)|(ipar == 5):
        x_parameter = np.linspace(Date_Num_parameter.min(), Date_Num_parameter.max(), 100)
        y1_parameter = np.linspace(depth_parameter.min(), depth_parameter.max(), 200)
        y2_parameter = np.linspace(dens_parameter.min(), dens_parameter.max(), 200)
        # I interpolate
        x_parameter_g, y_parameter_g = np.meshgrid(x_parameter, y1_parameter)
        parameter_interp_depth = griddata((Date_Num_parameter, depth_parameter), parameter_filtered,(x_parameter_g, y_parameter_g), method="nearest")
        x_parameter_g, y_parameter_g = np.meshgrid(x_parameter, y2_parameter)
        parameter_interp_dens = griddata((Date_Num_parameter,dens_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")

        sel200=y1_parameter<=200
        tmp=parameter_interp_depth[sel200,:]
        tmp=np.mean(tmp,axis=0)

        # Parameter in the top 200m
        fig = plt.figure(1, figsize=(12, 4))
        ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
        plt.plot(x_parameter,tmp)
        plt.ylabel(parameter_ylabel_list[ipar])
        plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
        plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
        plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
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
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        # I save
        plt.savefig('../Plots/an37/TimeSeriesPerDepthLayer/ZTimeSeries_%s_top200m_an37.pdf' % parameter_shortname_list[ipar],dpi=200)
        plt.close()

        # Parameter in the mixed layer
        mld_int = np.interp(x_parameter, Date_Num, mld)
        parameter_mld=np.zeros((mld_int.size,))
        i=0
        for i in range(0,mld_int.size):
            tmp=parameter_interp_depth[:,i]
            sel_mld=y1_parameter<=(mld_int[i]-20)
            parameter_mld[i]=np.mean(tmp[sel_mld])

        if ipar==1:
            temp_mld = parameter_mld.copy()
        elif ipar==4:
            parameter_mld = parameter_mld * mld_int
            chl_mld = parameter_mld.copy()
        elif ipar==5:
            doxy_mld = parameter_mld.copy()

        fig = plt.figure(1, figsize=(12, 4))
        ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
        plt.plot(x_parameter,parameter_mld)
        if ipar == 4:
            plt.ylabel('Integrated %s' % parameter_ylabel_list[ipar])
        else:
            plt.ylabel(parameter_ylabel_list[ipar])
        plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
        plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
        plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
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
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        # I save
        plt.savefig('../Plots/an37/ZTimeSeries_%s_ML_an37.pdf' % parameter_shortname_list[ipar],dpi=200)
        plt.close()

        if ipar ==4: continue
        # Parameter below mixed layer (1026.35-1027.24 kg/m3 isopycnals)
        parameter_mld_102724 = np.zeros((mld_int.size,))
        dens0_list=np.array([1026.35,1026.4,1026.8])
        dens0=dens0_list[0]
        dens1 = 1027.24  # 600 m
        for dens0 in dens0_list:
            i=0
            for i in range(0,mld_int.size):
                tmp=parameter_interp_dens[:,i]
                sel_tmp=(y2_parameter>=dens0)&(y2_parameter<dens1)
                parameter_mld_102724[i]=np.mean(tmp[sel_tmp])

            parameter_mld_102724 = savgol_filter(parameter_mld_102724, 5, 1)
            if ipar==1:
                temp_mld_102724 = parameter_mld_102724.copy()
            elif ipar==5:
                doxy_mld_102724 = parameter_mld_102724.copy()

            fig = plt.figure(1, figsize=(12, 4))
            ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
            plt.plot(x_parameter,parameter_mld_102724)
            plt.ylabel(parameter_ylabel_list[ipar])
            plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
            plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
            plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
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
            # I add the grid
            plt.grid(color='k', linestyle='dashed', linewidth=0.5)
            # I save
            plt.savefig('../Plots/an37/ZTimeSeries_%s_z%0.2f_%0.2fkgm3_an37.pdf' % (parameter_shortname_list[ipar],dens0,dens1),dpi=200)
            plt.close()


        # Parameter in the 200—600m layer
        sel200_600=(y1_parameter>200)&(y1_parameter<=600)
        tmp=parameter_interp_depth[sel200_600,:]
        tmp=np.mean(tmp,axis=0)
        if ipar==1:
            temp200_600m = tmp.copy()
        elif ipar==5:
            doxy200_600m = tmp.copy()

        fig = plt.figure(1, figsize=(12, 4))
        ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
        plt.plot(x_parameter,tmp)
        plt.ylabel(parameter_ylabel_list[ipar])
        plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
        plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
        plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
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
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        # I save
        plt.savefig('../Plots/an37/TimeSeriesPerDepthLayer/ZTimeSeries_%s_z200_600m_an37.pdf' % parameter_shortname_list[ipar],dpi=200)
        plt.close()




O2sat_interp = O2sat_func(np.reshape(psal_interp,(psal_interp.size,1)),np.reshape(cons_temp_interp_depth,(cons_temp_interp_depth.size,1)))
O2sat_interp = O2sat_interp.reshape(cons_temp_interp_depth.shape)

AOU_interp = O2sat_interp - doxy_interp
parameter_interp_depth=AOU_interp.copy()

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
plt.gca().invert_yaxis()
plt.ylim(ax.get_ylim()[1],ax.get_ylim()[0])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
# draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.set_ylabel('AOU ($\mu$mol/kg)', fontsize=18)
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
plt.savefig('../Plots/an37/TimeSeries_AOU_vs_01depth_an37.pdf',dpi=200)
plt.close()


# Parameter in the mixed layer: I approximate the AOU at 20m as representative of the ML, but should be improved in future
mld_int = np.interp(x_parameter, Date_Num, mld)
AOU_mld = np.zeros((mld_int.size,))
i = 0
for i in range(0, mld_int.size):
    tmp = AOU_interp[:, i]
    sel_mld = y1_parameter <= (mld_int[i]-20)
    AOU_mld[i] = np.mean(tmp[sel_mld])

# AOU_mld=AOU_interp[1,:]
fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(x_parameter,AOU_interp[1,:])
plt.ylabel('AOU ($\mu$mol/kg)')
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
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
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# I save
plt.savefig('../Plots/an37/ZTimeSeries_AOU_ML_an37.pdf' ,dpi=200)
plt.close()

#Mean chl in the mixed layer (and not as integrated chl a as plotted above)
width, height = 0.8, 0.7
mean_chl_mld = chl_mld/mld_int

fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.35, width, height-0.15])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(x_parameter,mean_chl_mld)
plt.ylabel('Chlorophyll-a (mg/m$^3$)')
plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
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
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# I save
plt.savefig('../Plots/an37/ZTimeSeries_chla_mean_ML_an37.pdf',dpi=200)
plt.close()


#######################################################################
# I save the mixed layer depth and BVfreq values for the latex document
#######################################################################
from write_latex_data import write_latex_data
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

Date_Vec_matlab=np.zeros((Date_Num.size,6))
for i in range(0,Date_Num.size):
    Date_Vec_matlab[i,:] = matlab_datevec(Date_Num[i]+matlab_datenum(1950,1,1))
Date_Vec_matlab = Date_Vec_matlab.astype(int)

filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
argument = 'mld_beginning'
arg_value=np.mean(mld[0:6])
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'mld_maximum'
arg_value=mld.max()
write_latex_data(filename,argument,'%d' % arg_value)
arg_value=Date_Vec_matlab[np.where(mld == mld.max()),:]
write_latex_data(filename,'mld_maximum_date','%d August' % arg_value[0][0][2])
argument = 'mld_end'
arg_value=mld[-1]
write_latex_data(filename,argument,'%d' % arg_value)

Date_Vec_parameter=np.zeros((x_parameter.size,6))
for i in range(0,x_parameter.size):
    Date_Vec_parameter[i,:] = matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1))
Date_Vec_parameter = Date_Vec_parameter.astype(int)

argument = 'bvf_75percentile_MLD'
arg_value=np.percentile(bvf_interp[0:5,:],75)*10**6
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'bvf_75percentile_exponent_MLD'
arg_value=-6
write_latex_data(filename,argument,'%d' % arg_value)

argument = 'bvf_102m_0413to0607'
arg_value=np.mean(bvf_interp[5,0:34])*10**4
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'bvf_102m_0413to0607_exponent'
arg_value=-4
write_latex_data(filename,argument,'%d' % arg_value)

argument = 'bvf_200to600m'
arg_value=np.mean(bvf_interp[10:-1,:])*10**6
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'bvf_200to600m_exponent'
arg_value=-6
write_latex_data(filename,argument,'%d' % arg_value)

argument = 'temp_ML_beginning'
arg_value=temp_mld[0]
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'temp_ML_end'
arg_value=temp_mld[-1]
write_latex_data(filename,argument,'%0.1f' % arg_value)

argument = 'doxy_ML_beginning'
arg_value=doxy_mld[0]
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'doxy_ML_end'
arg_value=doxy_mld[-1]
write_latex_data(filename,argument,'%d' % arg_value)

argument = 'AOU_ML_0413to0625'
arg_value=np.mean(AOU_mld[0:44])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'AOU_ML_0625to0905'
arg_value=np.mean(AOU_mld[44:89])
write_latex_data(filename,argument,'%0.1f' % arg_value)

argument = 'temp200_600m_0413'
arg_value=temp200_600m[0]
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'temp200_600m_maximum'
arg_value=temp200_600m.max()
write_latex_data(filename,argument,'%0.1f' % arg_value)
arg_value=matlab_datevec( x_parameter[np.where(temp200_600m == temp200_600m.max())][0] + matlab_datenum(1950,1,1) )
write_latex_data(filename,'temp200_600m_maximum_date','%d August' % np.floor( arg_value[2] ) )
argument = 'temp200_600m_0923'
arg_value=temp200_600m[-1]
write_latex_data(filename,argument,'%0.1f' % arg_value)

i_start=4
argument = 'temp102680_102721_startvalue'
arg_value=temp_mld_102724[i_start]
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'temp102680_102721_startdate'
arg_value=matlab_datevec(x_parameter[i_start]+matlab_datenum(1950,1,1)).astype(int)
write_latex_data(filename,argument,'%d April' % arg_value[2])
argument = 'temp102680_102721_endvalue'
arg_value=temp_mld_102724[-1]
write_latex_data(filename,argument,'%0.2f' % arg_value)

#######################################################################
# I save the mean and integrated chl concentration values for the latex document
#######################################################################
argument = 'max_chl_value'
arg_value=np.max(mean_chl_mld)
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'max_chl_date'
arg_value=np.where(mean_chl_mld==np.max(mean_chl_mld))[0][0]
arg_value=matlab_datevec((x_parameter[arg_value]+matlab_datenum(1950,1,1)))[2]
write_latex_data(filename,argument,'%d May' % arg_value)
i=73;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'min_chl_value'
arg_value=mean_chl_mld[i]
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'min_chl_date'
arg_value=matlab_datevec((x_parameter[i]+matlab_datenum(1950,1,1)))[2]
write_latex_data(filename,argument,'%d August' % arg_value)
i=90;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'increase_chl_date'
arg_value=matlab_datevec((x_parameter[i]+matlab_datenum(1950,1,1)))[2]
write_latex_data(filename,argument,'%d September' % arg_value)


argument = 'integrated_chl_202104'
arg_value=np.round(chl_mld[0])
write_latex_data(filename,argument,'%d' % arg_value)

i=69;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))

argument = 'integrated_chl_20210805'
arg_value=np.round(chl_mld[69])
write_latex_data(filename,argument,'%d' % arg_value)

