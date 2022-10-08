import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime,timedelta
from time import gmtime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import seawater as sw
import gsw
import netCDF4 as nc
import pickle
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec


#######################################################################
#######################################################################
#######################################################################
# I load the bbp data
#######################################################################
#######################################################################
#######################################################################
filename='6903095_Sprof_old.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])

Date_Num_bbp=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
date_reference = datetime.strptime("1/1/1950", "%d/%m/%Y")
date_reference_datenum=date.toordinal(date_reference)+366

Date_Vec=np.zeros([Date_Num_bbp.size,6])
for i in range(0,Date_Num_bbp.size):
    date_time_obj = date_reference + timedelta(days=Date_Num_bbp[i]-matlab_datenum(1950,1,1))
    Date_Vec[i,0]=date_time_obj.year;Date_Vec[i,1]=date_time_obj.month;Date_Vec[i,2]=date_time_obj.day
    Date_Vec[i,3]=date_time_obj.hour;Date_Vec[i,4]=date_time_obj.minute;Date_Vec[i,5]=date_time_obj.second

Date_Vec=Date_Vec.astype(int)

#Standard variables
pres=np.array(ds.variables['PRES_ADJUSTED'])
temp=np.array(ds.variables['TEMP_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])

#BGC Variables
bbp700=np.array(ds.variables['BBP700_ADJUSTED'])
chl = np.array(ds.variables['CHLA_ADJUSTED'])

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
if np.sum(bbp700==99999)==bbp700.size:
    print('Taking non adjusted bbp700')
    bbp700 = np.array(ds.variables['BBP700'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
    psal_qc = np.array(ds.variables['PSAL_QC'])
if np.sum(chl == 99999) == chl.size:
    print('Taking non adjusted chlorophyll-a')
    chl = np.array(ds.variables['CHLA'])

#######################################################################
#I tranform the pressure to depth
#######################################################################
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth_bbp=np.ones(pres.shape)*99999
depth_bbp[mask_depth]=depth_tmp

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
dens_bbp = np.ones(temp.shape) * 99999
dens_bbp[mask_dens] = dens_tmp + 1000

#######################################################################
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from oceanpy import bbp700toPOC
from oceanpy import bbp700toPOC_Koestner
bbp_POC=bbp700.copy()*0+99999
bbp_POC_Koestner=bbp700.copy()*0+99999
i=0
for i in range(0,bbp700.shape[0]):
    bbp700tmp=bbp700[i,:]
    depth_tmp=depth_bbp[i,:]
    temp_tmp=temp[i,:]
    chl_tmp=chl[i,:]
    # I exclude nan values
    sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)&(chl_tmp!=99999)
    bbp700tmp=bbp700tmp[sel]
    depth_tmp=depth_tmp[sel]
    temp_tmp=temp_tmp[sel]
    chl_tmp=chl_tmp[sel]
    # I convert to small POC (sPOC) and I set to 0 values <0
    sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
    sPOC_tmp[sPOC_tmp<0]=0
    bbp_POC[i,sel]=sPOC_tmp
    sPOC_tmp = bbp700toPOC_Koestner(bbp700tmp, chl_tmp)
    sPOC_tmp[np.isnan(sPOC_tmp)]=0
    bbp_POC_Koestner[i,sel]=sPOC_tmp

#######################################################################
# I convert the bbp dates to float values (in seconds from 1970 1 1)
#######################################################################
Date_Num_bbp_calendar = Date_Num_bbp.copy()
for i in range(0, Date_Num_bbp_calendar.size):
    date_time_obj = datetime(Date_Vec[i, 0], Date_Vec[i, 1], Date_Vec[i, 2],
                             Date_Vec[i, 3], Date_Vec[i, 4], Date_Vec[i, 5])
    Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.utcfromtimestamp(Date_Num[i])

#######################################################################
#######################################################################
#######################################################################
# I load the MiP MaP data
#######################################################################
#######################################################################
#######################################################################

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
lon=np.array(data['Longitude'][sel_filename])
lat=np.array(data['Latitude'][sel_filename])
Date_Time=np.array(data['Date_Time'][sel_filename])
pressure=np.array(data['Pressure [dbar]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
MiP_abund=np.array(data['MiP_abun'][sel_filename])
MaP_abund=np.array(data['MaP_abun'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
Date_Num_matlab=np.squeeze( np.zeros((Flux.size,1)) )
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    Date_Num_matlab[i] = matlab_datenum(date_time_obj.year,date_time_obj.month,date_time_obj.day,date_time_obj.hour,date_time_obj.minute,date_time_obj.second)
    #datetime.utcfromtimestamp(Date_Num[i])


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
Date_Num_Eddy_calendar=[]
i=Date_Num_Eddy[0]
for i in Date_Num_Eddy:
    daytmp = datetime.fromordinal(int(i)-366)
    DateTime_Eddy=np.append(DateTime_Eddy,daytmp)
    Date_Num_Eddy_calendar = np.append( Date_Num_Eddy_calendar,calendar.timegm(daytmp.timetuple()) )

#######################################################################
# I select the eddy position calculated from satellite in the days of the BGC Argo float profiles
#######################################################################
list_dates=np.sort(np.unique(Date_Num))

sel_Satel2Float=np.squeeze(np.zeros((1,list_dates.size)))
i=0
for i in range(0,list_dates.size):
    datetmp = gmtime(list_dates[i])
    datetmp = datetime(datetmp.tm_year,datetmp.tm_mon,datetmp.tm_mday)
    idx=np.array(np.where(datetmp == DateTime_Eddy))
    if idx.size>1: print('error, i %d, ' % i, idx)
    if idx.size==0:
        # print('Warning: missing eddy contour and center for %d-%d-%d' % (datetmp.year,datetmp.month,datetmp.day))
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

#I identify the unique list of lon and lat defining the BGC argo profile positions
lon_float=np.array([]);lat_float=np.array([])
i=0
for i in range(0, list_dates.size):
    sel = Date_Num == list_dates[i]
    lon_float = np.append(lon_float,np.unique(lon[sel]))
    lat_float = np.append(lat_float,np.unique(lat[sel]))

dist=np.arccos(np.sin(lat_float*np.pi/180)*np.sin(latEddy_4float*np.pi/180)+np.cos(lat_float*np.pi/180)*np.cos(latEddy_4float*np.pi/180)*np.cos((lonEddy_4float-lon_float)*np.pi/180))*180/np.pi
dist_km=dist*111


#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
# Before that, I save the bbpPOC and the depth so that it is easier to include it in other scripts
#######################################################################
sel_insideEddy = dist_km <= radius_Vmax_4float

dictionary_data = {"bbp_POC": bbp_POC,"bbp_POC_Koestner": bbp_POC_Koestner,"sel_insideEddy": sel_insideEddy,"Date_Num_bbp": Date_Num_bbp,
                   "Date_Vec_bbp": Date_Vec,"depth_bbp": depth_bbp,"Date_Num_bbp_calendar" : Date_Num_bbp_calendar,"temperature" : temp}
a_file = open("%s/an18/data_an18.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()

list_dates=list_dates[sel_insideEddy]
# lon=lon[sel_insideEddy]
# lat=lat[sel_insideEddy]
Date_Num_bbp=Date_Num_bbp[sel_insideEddy]
Date_Num_bbp_calendar=Date_Num_bbp_calendar[sel_insideEddy]
# Date_Time=Date_Time[sel_insideEddy]
depth_bbp=depth_bbp[sel_insideEddy]
temp=temp[sel_insideEddy]
# pressure=pressure[sel_insideEddy]
# Flux=Flux[sel_insideEddy]
# MiP_abund=MiP_abund[sel_insideEddy]
# MiP_POC=MiP_POC[sel_insideEddy]
# MaP_abund=MaP_abund[sel_insideEddy]
# MaP_POC=MaP_POC[sel_insideEddy]
bbp_POC=bbp_POC[sel_insideEddy,:]

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from oceanpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,temp.shape[0]):
    depth_tmp=depth_bbp[i,:]
    temp_tmp=temp[i,:]
    # I exclude nan values
    sel_non_nan=(depth_tmp!=99999)&(temp_tmp!=99999)
    temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
    mld_tmp,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
    mld=np.append(mld,mld_tmp)

#######################################################################
# I start the loop on the different parameters I plot
#######################################################################
parameter=Flux
ipar=5
parameter_shortname_list=['Flux','MiP_abund','MaP_abund','MiP_POC','MaP_POC','bbpPOC']
parameter_panellabel_list=['b','g','h','g','h','f']
parameter_ylabel_list=['Flux (mgC $m^{-2}$ $d^{-1}$)','MiP abundance (# L$^{-1}$)','MaP abundance (# L$^{-1}$)'
    ,'MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)','$b_{bp}$POC (mgC $m^{-3}$)']
max_parameter_list=np.array([32,65,0.6,2.15,0.30,40])
MiP_POC_0_200=np.array([]);MiP_POC_200_600=np.array([])
MaP_POC_0_200=np.array([]);MaP_POC_200_600=np.array([])
bbp_POC_0_200=np.array([]);bbp_POC_200_600=np.array([])
max_depth_bbp=np.array([]);max_depth=np.array([])
day_start_eddy_merging = datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0: parameter=Flux.copy()
    elif ipar==1:   parameter=MiP_abund.copy()
    elif ipar == 2: parameter=MaP_abund.copy()
    elif ipar == 3: parameter=MiP_POC.copy()
    elif ipar == 4: parameter=MaP_POC.copy()
    elif ipar == 5: parameter=bbp_POC.copy()

    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([]);max_depth_parameter=np.array([]);
    if ipar == 5:
        i=0
        for i in range(0, bbp_POC.shape[0]):
            z=parameter[i,:];y=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
            max_depth = np.append(max_depth ,np.unique(np.sort(y))[-2]) #-2 to exclude 99999
            z[z>100] = 99999
            sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
            sel3 = z == 0
            if sum(sel2) > 0:
                z = savgol_filter(z, 5, 1)
                z[sel3] = 0
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, np.tile(x,sum(sel2)) ))
                depth_filtered = np.concatenate((depth_filtered, y2))
                max_depth_bbp = np.append(max_depth_bbp,y2.max())
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
                # if (ipar==5)&(i==36):   z[-1]=0 #With this line, I exclude a spike measured in one bbpPOC profile at 600 m which was making the bbpPOC integrated time series odd (in the sense that it had a anamalous spike corresponding to that profile)
                z=savgol_filter(z,5,1)
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
                depth_filtered = np.concatenate((depth_filtered, y2))
                max_depth_parameter = np.append(max_depth_parameter, y2.max())
                # sel_200 and sel_200_600 are used only for the POC integrated in time
                sel_0_200 = np.abs(y2) < 200
                sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
                if ipar==3: MiP_POC_0_200=np.append(MiP_POC_0_200,np.mean(z[sel_0_200]));MiP_POC_200_600=np.append(MiP_POC_200_600,np.mean(z[sel_200_600]))
                if ipar==4: MaP_POC_0_200=np.append(MaP_POC_0_200,np.mean(z[sel_0_200]));MaP_POC_200_600=np.append(MaP_POC_200_600,np.mean(z[sel_200_600]))
                # if ipar==5: bbp_POC_0_200=np.append(bbp_POC_0_200,np.mean(z[sel_0_200]));bbp_POC_200_600=np.append(bbp_POC_200_600,np.mean(z[sel_200_600]))

    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    sel_0_200 = (np.abs(y_filtered) >= 0) & (np.abs(y_filtered) < 200)
    sel_200_600 = (np.abs(y_filtered) >= 200) & (np.abs(y_filtered) < 600)
    if ipar==3: MiP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MiP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==4: MaP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MaP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==5: bbp_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);bbp_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)

    if ipar==3: MIP_POC_interp = parameter_interp.copy();max_depth_MIP=max_depth_parameter
    if ipar==4: MAP_POC_interp = parameter_interp.copy();max_depth_MAP=max_depth_parameter
    if ipar==5: bbp_POC_interp = parameter_interp.copy()

    ################################################################################################################
    ####### I plot
    ################################################################################################################

    # if ipar == 5: continue  # I do not plot the bbp as it is already plotted in an17 (with higher resolution data)

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w')
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
        xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90,fontsize=12)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an18/TimeSeries%02dA%s_an18.pdf' % (ipar,parameter_shortname_list[ipar]),dpi=200)
    plt.close()

    ####### I do the interpolate contourf up to the maximum depth
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),1000
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=6000,color='w')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=6000,color='w')
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
        xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90,fontsize=12)
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an18/TimeSeries%02dB%s_an18.pdf' % (ipar,parameter_shortname_list[ipar]),dpi=200)
    plt.close()



###########################################################################
###########################################################################
# I plot the integrated POC budget over time, for the layer 0-200m and 200-600m
###########################################################################
###########################################################################

POC_0_200=MiP_POC_0_200+MaP_POC_0_200+bbp_POC_0_200
POC_0_200_int=MiP_POC_0_200_int+MaP_POC_0_200_int+bbp_POC_0_200_int
POC_200_600=MiP_POC_200_600+MaP_POC_200_600+bbp_POC_200_600
POC_200_600_int=MiP_POC_200_600_int+MaP_POC_200_600_int+bbp_POC_200_600_int
POC_0_200[POC_0_200<0]=0
POC_0_200_int[POC_0_200_int<0]=0
POC_200_600[POC_200_600<0]=0
POC_200_600_int[POC_200_600_int<0]=0
# Parameters for the plot
width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(POC_0_200.min(),POC_200_600.min()*10),max(POC_0_200.max(),POC_200_600.max()*10)

#POC 0-200 vs 200-600
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
# plt.plot(list_dates,POC_0_200,'r',label='0-200 m')
# plt.scatter(list_dates,POC_0_200,c='r')
# plt.plot(list_dates,POC_200_600*10,'b',label='200-600 m')
# plt.scatter(list_dates,POC_200_600*10,c='b')
plt.plot(x_filtered,POC_0_200_int,'r',linewidth=3,label='0-200 m')
plt.plot(x_filtered,POC_200_600_int*10,'b',linewidth=3,label='200-600 m')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
ax.text(-0.075, 1.05, 'e', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an18/IntegratedPOC_vs_time_an18.pdf' ,dpi=200)
plt.close()



# 0-200 m layer
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,bbp_POC_0_200,'y',label='bbpPOC')
plt.scatter(list_dates,bbp_POC_0_200,c='y')
plt.plot(list_dates,MiP_POC_0_200,'c',label='MiP')
plt.scatter(list_dates,MiP_POC_0_200,c='c')
plt.plot(list_dates,MaP_POC_0_200,'g',label='MaP')
plt.scatter(list_dates,MaP_POC_0_200,c='g')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
plt.title('0-200 m layer', fontsize=18)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an18/0_200m_MiP_MaP_bbpPOC_vs_time_an18.pdf' ,dpi=200)
plt.close()




# 200-600 m layer
set_ylim_lower, set_ylim_upper = POC_200_600.min(),POC_200_600_int.max()
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
# plt.plot(list_dates,bbp_POC_200_600,'y',label='bbpPOC')
# plt.scatter(list_dates,bbp_POC_200_600,c='y')
# plt.plot(list_dates,MiP_POC_200_600,'c',label='MiP')
# plt.scatter(list_dates,MiP_POC_200_600,c='c')
# plt.plot(list_dates,MaP_POC_200_600,'g',label='MaP')
# plt.scatter(list_dates,MaP_POC_200_600,c='g')
plt.plot(x_filtered,bbp_POC_200_600_int,'y',linewidth=2,label='bbpPOC')
plt.plot(x_filtered,MiP_POC_200_600_int,'c',linewidth=2,label='MiP')
plt.plot(x_filtered,MaP_POC_200_600_int,'g',linewidth=2,label='MaP')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
plt.title('200-600 m layer', fontsize=18)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an18/200_600m_MiP_MaP_bbpPOC_vs_time_an18.pdf' ,dpi=200)
plt.close()


#######################################################################
# I save the MIP, MAP, and POC concentration values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum

datetime.utcfromtimestamp(x_filtered[44])
argument = 'POC_0_200m_0413to0625'
arg_value=np.mean(POC_0_200_int[0:45])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'POC_0_200m_maximum'
arg_value=POC_0_200_int.max()
write_latex_data(filename,argument,'%0.1f' % arg_value)
arg_value=datetime.utcfromtimestamp( x_filtered[np.where(POC_0_200_int == POC_0_200_int.max())][0] )
write_latex_data(filename,'POC_0_200m_maximum_date','%d August' %  arg_value.day  )
peak_POC_0_200m_date=arg_value.day

argument = 'POC_200_600m_0413'
arg_value=POC_200_600_int[0]
write_latex_data(filename,argument,'%0.1f' % arg_value)
datetime.utcfromtimestamp(x_filtered[44])
argument = 'POC_200_600m_0625'
arg_value=POC_200_600_int[44]
write_latex_data(filename,argument,'%0.1f' % arg_value)
POC_200_600_int[66:74]
arg_value = datetime.utcfromtimestamp(x_filtered[70]).day
argument = 'POC_200_600m_before_peak'
write_latex_data(filename,argument,'%d August' % arg_value)
datetime.utcfromtimestamp(x_filtered[73])
POC_200_600_int[71:78]
peak_POC_200_600m_date=datetime.utcfromtimestamp(x_filtered[73]).day
argument = 'POC_200_600m_08%02d' % peak_POC_200_600m_date
arg_value=POC_200_600_int[73]
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'POC_200_600m_08%02d_date' % peak_POC_200_600m_date
write_latex_data(filename,argument,'%d August' % peak_POC_200_600m_date)
argument = 'peak_0_200_200_600_mismatch'
arg_value = peak_POC_200_600m_date - peak_POC_0_200m_date
write_latex_data(filename,argument,'%d' % arg_value)

y_filtered[9]
datetime.utcfromtimestamp(x_filtered[22])
argument = 'Map_POC_2021AprMay'
arg_value=np.mean(MAP_POC_interp[0:9,0:23])
write_latex_data(filename,argument,'%0.2f' % arg_value)
datetime.utcfromtimestamp(x_filtered[66])
datetime.utcfromtimestamp(x_filtered[88])
y_filtered[16]
argument = 'Map_POC_2021FluxEvent'
arg_value=np.mean(MAP_POC_interp[0:16,66:89])
write_latex_data(filename,argument,'%0.2f' % arg_value)

argument = 'Mip_POC_2021AprMay'
arg_value=np.mean(MIP_POC_interp[0:9,0:23])
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Mip_POC_2021FluxEvent'
arg_value=np.mean(MIP_POC_interp[0:16,66:89])
write_latex_data(filename,argument,'%0.2f' % arg_value)
y_filtered[16]
argument = 'Mip_POC_0413_90to600m'
arg_value=np.mean(MIP_POC_interp[9:60,0])
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Mip_POC_2021FluxEvent_160to600m'
arg_value=np.mean(MIP_POC_interp[16:60,66:89])
write_latex_data(filename,argument,'%0.2f' % arg_value)

argument = 'Mip_POC_0413_200to600m'
arg_value=np.mean(MiP_POC_200_600_int[0])
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Mip_POC_200to600m_min'
arg_value=MiP_POC_200_600_int.min()
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'Mip_POC_200to600m_min_date'
arg_value = datetime.utcfromtimestamp(x_filtered[ np.where(MiP_POC_200_600_int == MiP_POC_200_600_int.min() )[0][0] ]).day
write_latex_data(filename,argument,'%d July' % arg_value)

#######################################################################
# I calculate the bbp and MaP POC in and outside the mld, I plot and save it for the latex document
#######################################################################
mld_int = np.interp(x_filtered,Date_Num_bbp_calendar,mld)
bbp_POC_mld=np.array([])
bbp_POC_out_mld=np.array([])
MaP_POC_mld=np.array([])
MaP_POC_out_mld=np.array([])
MiP_POC_mld=np.array([])
MiP_POC_out_mld=np.array([])
i=0
for i in range(0,bbp_POC_interp.shape[1]):
    sel_tmp = y_filtered<=mld_int[i]
    bbp_POC_mld = np.append(bbp_POC_mld, np.mean( bbp_POC_interp[sel_tmp,i] ))
    bbp_POC_out_mld = np.append(bbp_POC_out_mld, np.mean( bbp_POC_interp[~sel_tmp,i] ))
    MaP_POC_mld = np.append(MaP_POC_mld, np.mean( MAP_POC_interp[sel_tmp,i] ))
    MaP_POC_out_mld = np.append(MaP_POC_out_mld, np.mean( MAP_POC_interp[~sel_tmp,i] ))
    MiP_POC_mld = np.append(MiP_POC_mld, np.mean( MIP_POC_interp[sel_tmp,i] ))
    MiP_POC_out_mld = np.append(MiP_POC_out_mld, np.mean( MIP_POC_interp[~sel_tmp,i] ))

fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.4, width, height])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(x_filtered,bbp_POC_mld,label='$b_{bp}$ POC in ML')
plt.plot(x_filtered,bbp_POC_out_mld*10,label='$b_{bp}$ POC out ML $\cdot$10')
plt.ylabel('bbp POC (mgC/m$^3$)')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.legend()
plt.ylim([0,bbp_POC_mld.max()*1.1])
# I save
plt.savefig('../Plots/an18/ZTimeSeries_bbp_ML_an18.pdf' ,dpi=200)
plt.close()

fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.4, width, height])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(x_filtered,MaP_POC_mld,label='MaP POC in ML')
plt.plot(x_filtered,MaP_POC_out_mld,label='MaP POC out ML')
plt.ylabel('MaP POC (mgC/m$^3$)')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.legend()
plt.ylim([0,MaP_POC_mld.max()*1.1])
# I save
plt.savefig('../Plots/an18/ZTimeSeries_MaP_ML_an18.pdf' ,dpi=200)
plt.close()

fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_axes([0.12, 0.4, width, height])# ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
plt.plot(x_filtered,MiP_POC_mld,label='MiP POC in ML')
plt.plot(x_filtered,MiP_POC_out_mld,label='MiP POC out ML')
plt.ylabel('MiP POC (mgC/m$^3$)')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.legend()
plt.ylim([0,MiP_POC_mld.max()*1.1])
# I save
plt.savefig('../Plots/an18/ZTimeSeries_MiP_ML_an18.pdf' ,dpi=200)
plt.close()

argument = 'bbp_ML_0413'
arg_value=bbp_POC_mld[0]
write_latex_data(filename,argument,'%0.1f' % arg_value)
datetime.utcfromtimestamp(x_filtered[44])
argument = 'bbp_ML_0625'
arg_value=np.mean(bbp_POC_mld[44])
write_latex_data(filename,argument,'%0.1f' % arg_value)
datetime.utcfromtimestamp(x_filtered[70])
bbp_POC_mld[65:75]
argument = 'bbp_ML_0807'
arg_value=np.mean(bbp_POC_mld[70])
write_latex_data(filename,argument,'%0.1f' % arg_value)
argument = 'bbp_outML_0413'
arg_value=bbp_POC_out_mld[0]
write_latex_data(filename,argument,'%0.1f' % arg_value)
datetime.utcfromtimestamp(x_filtered[48])
argument = 'bbp_outML_0701to0923'
arg_value=np.mean(bbp_POC_out_mld[48:])
write_latex_data(filename,argument,'%0.1f' % arg_value)

argument = 'MaP_outML_0413to0923'
arg_value=np.mean(MaP_POC_out_mld)
write_latex_data(filename,argument,'%0.02f' % arg_value)
argument = 'MaP_outML_0413to0923_std'
arg_value=np.std(MaP_POC_out_mld)
write_latex_data(filename,argument,'%0.02f' % arg_value)

i=59;print(datetime.utcfromtimestamp(x_filtered[i]).strftime('%d %B'))
argument = 'MiP_ML_0413to0719'
arg_value=np.mean(MiP_POC_mld[0:59])
write_latex_data(filename,argument,'%0.02f' % arg_value)
i=71;print(datetime.utcfromtimestamp(x_filtered[i]).strftime('%d %B'))
arg_value=np.mean(MiP_POC_mld[i])
argument = 'MiP_ML_0808'
write_latex_data(filename,argument,'%0.02f' % arg_value)

#I calculate also the difference of the integrated POC 200—600 m between the 13 April and the 20 June 2021 (it is a statistic for the main paper)
datetime.utcfromtimestamp(x_filtered[41])
argument = 'Integrated_POC_0620_0413difference'
arg_value = (POC_200_600_int[0]-POC_200_600_int[41])*400
# write_latex_data(filename,argument,'%d' % arg_value) #I don't save it anymore cos I use an53

a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
data_an18 = pickle.load(a_file)
bbp_POC = data_an18['bbp_POC']
bbp_POC_Koestner = data_an18['bbp_POC_Koestner']
sel_insideEddy = data_an18['sel_insideEddy']
Date_Num_bbp = data_an18['Date_Num_bbp']
Date_Vec_bbp = data_an18['Date_Vec_bbp']
Date_Num_bbp_calendar = data_an18['Date_Num_bbp_calendar']
depth_bbp = data_an18['depth_bbp']
a_file.close()

dictionary_data = {"bbp_POC": bbp_POC,"bbp_POC_Koestner": bbp_POC_Koestner,"sel_insideEddy": sel_insideEddy,"Date_Num_bbp": Date_Num_bbp,
                   "Date_Vec_bbp": Date_Vec,"depth_bbp": depth_bbp,"Date_Num_bbp_calendar": Date_Num_bbp_calendar,
                   "Integrated_POC_0620_0413difference": arg_value, "dens_bbp": dens_bbp}
a_file = open("%s/an18/data_an18.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()

#I plot the maximum depth
width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = 0,max_depth.max()
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,max_depth,'.b-',label='abs.')
plt.plot(list_dates,max_depth_bbp,'.r-',label='bbp')
plt.plot(list_dates,max_depth_MIP,'.m-',label='MiP')
plt.plot(list_dates,max_depth_MAP,'.c-',label='MaP')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=6000, color='k')
plt.vlines(day_end_eddy_merging, ymin=0, ymax=6000, color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('depth (m)', fontsize=15)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an18/Max_depth_vs_profile_an18.pdf' ,dpi=200)
plt.close()

#I compare interpolated and not interpolated profiles one by one
fs=10
width, height = 0.76, 0.7
nprofiles_to_plot=12 # I plot only the first 12 profiles
ipar=5
for ipar in range(3,parameter_ylabel_list.__len__()):
    if ipar == 3: parameter=MiP_POC.copy()
    elif ipar == 4: parameter=MaP_POC.copy()
    elif ipar == 5: parameter=bbp_POC.copy()

    if ipar == 5:
        i=0
        for i in range(0, nprofiles_to_plot):
            z=parameter[i,:];y=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
            max_depth = np.append(max_depth ,np.unique(np.sort(y))[-2]) #-2 to exclude 99999
            z[z>100] = 99999
            sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
            sel3 = z == 0
            if sum(sel2) > 0:
                z = savgol_filter(z, 5, 1)
                z[sel3] = 0
                # The bbp is interpolated at dates which do not correspond with the exact dates of the profiles: for this
                # reason, I take the interpolated profile whose date is the closest one to the real profile considered
                # in the loop
                idx=abs(x_filtered-x)
                idx = np.where(idx==idx.min())[0][0]
                datetmp = datetime.utcfromtimestamp(x)
                fig = plt.figure(1, figsize=(3.5, 3.5))
                ax = fig.add_axes([0.2, 0.15, width, height])
                plt.plot(bbp_POC_interp[sel_200_600,idx],y_filtered[sel_200_600],'.b-',label='interp.')
                plt.plot(z[seltmp],y2[seltmp],'.g-',label='raw')
                plt.xlim([-0.2,10])
                plt.ylim([195,600])
                plt.gca().invert_yaxis()
                plt.xlabel(parameter_ylabel_list[ipar], fontsize=fs)
                plt.ylabel('Depth (m)', fontsize=fs)
                plt.title('Profile %d, %s\n200-600m POC, raw:%0.2f int:%0.2f' % (int(i)+1,datetmp.__str__(),bbp_POC_200_600[i],bbp_POC_200_600_int[idx]))
                plt.legend(fontsize=fs)
                plt.grid(color='k', linestyle='dashed', linewidth=0.5)
                plt.savefig('../Plots/an18/Interpolated_vs_not/%s_profile%02d_an18.pdf' % (parameter_shortname_list[ipar],int(i)+1), dpi=200)
                plt.close()



    else:
        if ipar==3: tmp = MIP_POC_interp.copy();tmp2=MiP_POC_200_600.copy();tmp3=MiP_POC_200_600_int.copy()
        if ipar==4: tmp = MAP_POC_interp.copy();tmp2=MaP_POC_200_600.copy();tmp3=MaP_POC_200_600_int.copy();max_parameter_list[ipar]=0.8
        # I filter the prophiles
        i=0
        for i in range(0,nprofiles_to_plot):
            sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
            z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
            if sum(sel2)>0:
                z=savgol_filter(z,5,1)
                # The Mip and Map are interpolated at dates which do not correspond with the exact dates of the profiles: for this
                # reason, I take the interpolated profile whose date is the closest one to the real profile considered
                # in the loop
                idx=abs(x_filtered-list_dates[i])
                idx = np.where(idx==idx.min())[0][0]
                datetmp = datetime.utcfromtimestamp(list_dates[i])
                fig = plt.figure(1, figsize=(3.5, 3.5))
                ax = fig.add_axes([0.2, 0.15, width, height])
                plt.plot(tmp[:,idx],y_filtered,'.b-',label='interp.')
                plt.plot(z,y2,'.g-',label='raw.')
                plt.xlim([0,max_parameter_list[ipar]])
                plt.ylim([195,600])
                plt.gca().invert_yaxis()
                plt.xlabel(parameter_ylabel_list[ipar], fontsize=fs)
                plt.ylabel('Depth (m)', fontsize=fs)
                plt.title('Profile %d, %s\n200-600m POC, raw:%0.2f int:%0.2f' % (int(i)+1,datetmp.__str__(),tmp2[i],tmp3[idx]))
                plt.legend(fontsize=fs)
                plt.grid(color='k', linestyle='dashed', linewidth=0.5)
                plt.savefig('../Plots/an18/Interpolated_vs_not/%s_profile%02d_an18.pdf' % (parameter_shortname_list[ipar],int(i)+1), dpi=200)
                plt.close()



a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
data_an18 = pickle.load(a_file)
bbp_POC = data_an18['bbp_POC']
sel_insideEddy = data_an18['sel_insideEddy']
Date_Num_bbp = data_an18['Date_Num_bbp']
Date_Vec_bbp = data_an18['Date_Vec_bbp']
depth_bbp = data_an18['depth_bbp']
Date_Num_bbp_calendar = data_an18['Date_Num_bbp_calendar']
Integrated_POC_0620_0413difference = data_an18['Integrated_POC_0620_0413difference']
dens_bbp = data_an18['dens_bbp']
a_file.close()

dictionary_data = {"bbp_POC": bbp_POC,"bbp_POC_Koestner": bbp_POC_Koestner,"sel_insideEddy": sel_insideEddy,"Date_Num_bbp": Date_Num_bbp,
                   "Date_Vec_bbp": Date_Vec,"depth_bbp": depth_bbp,"Date_Num_bbp_calendar": Date_Num_bbp_calendar,
                   "Integrated_POC_0620_0413difference": Integrated_POC_0620_0413difference, "dens_bbp": dens_bbp}
a_file = open("%s/an18/data_an18.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()
