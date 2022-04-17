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
filename='6903095_Sprof.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])

Date_Num_bbp=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
date_reference = datetime.strptime("1/1/1950", "%d/%m/%Y")
date_reference_datenum=date.toordinal(date_reference)+366

Date_Vec=np.zeros([Date_Num_bbp.size,6])
for i in range(0,Date_Num_bbp.size):
    date_time_obj = date_reference + timedelta(days=Date_Num_bbp[i])
    Date_Vec[i,0]=date_time_obj.year;Date_Vec[i,1]=date_time_obj.month;Date_Vec[i,2]=date_time_obj.day
    Date_Vec[i,3]=date_time_obj.hour;Date_Vec[i,4]=date_time_obj.minute;Date_Vec[i,5]=date_time_obj.second

Date_Vec=Date_Vec.astype(int)

#Standard variables
pres=np.array(ds.variables['PRES_ADJUSTED'])
temp=np.array(ds.variables['TEMP_ADJUSTED'])

#BGC Variables
bbp700=np.array(ds.variables['BBP700_ADJUSTED'])

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
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from paruvpy import bbp700toPOC
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

dictionary_data = {"bbp_POC": bbp_POC,"sel_insideEddy": sel_insideEddy,"Date_Num_bbp": Date_Num_bbp,
                   "Date_Vec_bbp": Date_Vec,"depth_bbp": depth_bbp,"sel_insideEddy": sel_insideEddy}
a_file = open("%s/an18/data_an18.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()


list_dates=list_dates[sel_insideEddy]

# lon=lon[sel_insideEddy]
# lat=lat[sel_insideEddy]
Date_Num_bbp=Date_Num_bbp[sel_insideEddy]
# Date_Time=Date_Time[sel_insideEddy]
depth_bbp=depth_bbp[sel_insideEddy]
# pressure=pressure[sel_insideEddy]
# Flux=Flux[sel_insideEddy]
# MiP_abund=MiP_abund[sel_insideEddy]
# MiP_POC=MiP_POC[sel_insideEddy]
# MaP_abund=MaP_abund[sel_insideEddy]
# MaP_POC=MaP_POC[sel_insideEddy]
bbp_POC=bbp_POC[sel_insideEddy,:]

#######################################################################
# I start the loop on the different parameters I plot
#######################################################################
parameter=Flux
ipar=1
parameter_shortname_list=['Flux','MiP_abund','MaP_abund','MiP_POC','MaP_POC']
parameter_panellabel_list=['b','g','h','g','h']
parameter_ylabel_list=['Flux (mgC $m^{-2}$ $d^{-1}$)','MiP abundance (# L$^{-1}$)','MaP abundance (# L$^{-1}$)'
    ,'MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)']
max_parameter_list=np.array([32,65,0.6,2.15,0.30])
MiP_POC_0_200=np.array([]);MiP_POC_200_600=np.array([])
MaP_POC_0_200=np.array([]);MaP_POC_200_600=np.array([])
bbp_POC_0_200=np.array([]);bbp_POC_200_600=np.array([])
for ipar in range(0,parameter_ylabel_list.__len__()+1):
    if ipar==0: parameter=Flux.copy()
    elif ipar==1:   parameter=MiP_abund.copy()
    elif ipar == 2: parameter=MaP_abund.copy()
    elif ipar == 3: parameter=MiP_POC.copy()
    elif ipar == 4: parameter=MaP_POC.copy()
    elif ipar == 5: parameter=bbp_POC.copy()

    if ipar == 5:
        i=0
        for i in range(0, bbp_POC.shape[0]):
            z=parameter[i,:];x=Date_Num_bbp[i];y=depth_bbp[i,:]
            z[z>100] = 99999
            sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
            if sum(sel2) > 0:
                z = savgol_filter(z, 5, 1)
                # sel_200 and sel_200_600 are used only for the POC integrated in time
                sel_0_200 = np.abs(y2) < 200
                sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
                bbp_POC_0_200=np.append(bbp_POC_0_200,np.mean(z[sel_0_200]));bbp_POC_200_600=np.append(bbp_POC_200_600,np.mean(z[sel_200_600]))
        continue

    # I filter the flux prophiles
    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
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
            # sel_200 and sel_200_600 are used only for the POC integrated in time
            sel_0_200 = np.abs(y2) < 200
            sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
            if ipar==3: MiP_POC_0_200=np.append(MiP_POC_0_200,np.mean(z[sel_0_200]));MiP_POC_200_600=np.append(MiP_POC_200_600,np.mean(z[sel_200_600]))
            if ipar==4: MaP_POC_0_200=np.append(MaP_POC_0_200,np.mean(z[sel_0_200]));MaP_POC_200_600=np.append(MaP_POC_200_600,np.mean(z[sel_200_600]))
            # if ipar==5: bbp_POC_0_200=np.append(bbp_POC_0_200,np.mean(z[sel_0_200]));bbp_POC_200_600=np.append(bbp_POC_200_600,np.mean(z[sel_200_600]))

    if ipar==5: continue # I do not plot the bbp as it is already plotted in an17 (with higher resolution data)
    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    ################################################################################################################
    ####### I plot
    ################################################################################################################
    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
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
POC_200_600=MiP_POC_200_600+MaP_POC_200_600+bbp_POC_200_600
POC_0_200[POC_0_200<0]=0
POC_200_600[POC_200_600<0]=0
# Parameters for the plot
width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(POC_0_200.min(),POC_200_600.min()*10),max(POC_0_200.max(),POC_200_600.max()*10)

#POC 0-200 vs 200-600
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,POC_0_200,'r',label='0-200 m')
plt.scatter(list_dates,POC_0_200,c='r')
plt.plot(list_dates,POC_200_600*10,'b',label='200-600 m')
plt.scatter(list_dates,POC_200_600*10,c='b')
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
set_ylim_lower, set_ylim_upper = POC_200_600.min(),POC_200_600.max()
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,bbp_POC_200_600,'y',label='bbpPOC')
plt.scatter(list_dates,bbp_POC_200_600,c='y')
plt.plot(list_dates,MiP_POC_200_600,'c',label='MiP')
plt.scatter(list_dates,MiP_POC_200_600,c='c')
plt.plot(list_dates,MaP_POC_200_600,'g',label='MaP')
plt.scatter(list_dates,MaP_POC_200_600,c='g')
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


