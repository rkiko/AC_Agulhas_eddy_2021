import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime
from time import gmtime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import seawater as sw
import pickle
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datevec import matlab_datevec

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
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
MiP_abund=np.array(data['MiP_abun'][sel_filename])
MaP_abund=np.array(data['MaP_abun'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])
bbp_POC=np.array(data['bbp POC [mgC/m3]'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
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
#######################################################################
sel_insideEddy = dist_km <= radius_Vmax_4float

list_dates=list_dates[sel_insideEddy]

#######################################################################
# I load the euphotic layer depth
#######################################################################
a_file = open("%s/an36/data_an36.pkl" % storedir, "rb")
data_an36 = pickle.load(a_file)
zeu_float=data_an36['zeu_float']
zeu_datenum=data_an36['zeu_datenum']
zeu_matlab_datenum=data_an36['zeu_matlab_datenum']

i=0
for i in range(0,zeu_datenum.size):
    if np.isnan(zeu_matlab_datenum[i]): continue
    datetmp = matlab_datevec(zeu_matlab_datenum[i])
    date_time_obj = datetime.strptime('%d-%02d-%02dT%02d:%02d:%02d' % (datetmp[0],datetmp[1],datetmp[2],
                    datetmp[3],datetmp[4],datetmp[5]), '%Y-%m-%dT%H:%M:%S')
    zeu_datenum[i] = calendar.timegm(date_time_obj.timetuple())

zeu_datenum = zeu_datenum[~np.isnan(zeu_float)]
zeu_float = zeu_float[~np.isnan(zeu_float)]

#######################################################################
# I start the loop on the different parameters I plot
#######################################################################
parameter=Flux
ipar=1
parameter_shortname_list=['Flux','Flux_eta_b','Flux_extended','Flux_extended_eta_b','MiP_abund','MaP_abund','MiP_POC','MiP_POC_extended','MaP_POC']
parameter_panellabel_list=['b','b','b','b','f','h','f','f','h']
parameter_ylabel_list=['Flux (mgC $m^{-2}$ $d^{-1}$)','Flux (mgC $m^{-2}$ $d^{-1}$)','Flux (mgC $m^{-2}$ $d^{-1}$)',
                       'Flux (mgC $m^{-2}$ $d^{-1}$)','MiP abundance (# L$^{-1}$)','MaP abundance (# L$^{-1}$)'
                        ,'MiP (mgC $m^{-3}$)','MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)']
max_parameter_list=np.array([32,16,100,50,65,0.6,2.15,25.15,0.30])
MiP_POC_0_200=np.array([]);MiP_POC_200_600=np.array([])
MiP_POC_extended_0_200=np.array([]);MiP_POC_extended_200_600=np.array([])
MaP_POC_0_200=np.array([]);MaP_POC_200_600=np.array([])
bbp_POC_0_200=np.array([]);bbp_POC_200_600=np.array([])
for ipar in range(0,parameter_ylabel_list.__len__()+1):
    if ipar==0: parameter=Flux.copy()
    if ipar==1: parameter=Flux_eta_b.copy()
    if ipar==2: parameter=Flux_extended.copy()
    if ipar==3: parameter=Flux_extended_eta_b.copy()
    elif ipar==4:   parameter=MiP_abund.copy()
    elif ipar == 5: parameter=MaP_abund.copy()
    elif ipar == 6: parameter=MiP_POC.copy()
    elif ipar == 7: parameter=MiP_POC_extended.copy()
    elif ipar == 8: parameter=MaP_POC.copy()
    elif ipar == 9: parameter=bbp_POC.copy()

    # I filter the flux prophiles
    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
        z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2)>0:
            if (ipar==9)&(i==36):   z[-1]=0 #With this line, I exclude a spike measured in one bbpPOC profile at 600 m which was making the bbpPOC integrated time series odd (in the sense that it had a anamalous spike corresponding to that profile)
            z=savgol_filter(z,5,1)
            parameter_filtered = np.concatenate((parameter_filtered, z))
            Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
            depth_filtered = np.concatenate((depth_filtered, y2))
            # sel_200 and sel_200_600 are used only for the POC integrated in time
            sel_0_200 = np.abs(y2) < 200
            sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
            if ipar==6: MiP_POC_0_200=np.append(MiP_POC_0_200,np.mean(z[sel_0_200]));MiP_POC_200_600=np.append(MiP_POC_200_600,np.mean(z[sel_200_600]))
            if ipar==7: MiP_POC_extended_0_200=np.append(MiP_POC_extended_0_200,np.mean(z[sel_0_200]));MiP_POC_extended_200_600=np.append(MiP_POC_extended_200_600,np.mean(z[sel_200_600]))
            if ipar==8: MaP_POC_0_200=np.append(MaP_POC_0_200,np.mean(z[sel_0_200]));MaP_POC_200_600=np.append(MaP_POC_200_600,np.mean(z[sel_200_600]))
            if ipar==9: bbp_POC_0_200=np.append(bbp_POC_0_200,np.mean(z[sel_0_200]));bbp_POC_200_600=np.append(bbp_POC_200_600,np.mean(z[sel_200_600]))

    if ipar==9: continue # I do not plot the bbp as it is already plotted in an17 (with higher resolution data)
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
    plt.plot(zeu_datenum, zeu_float, 'w')
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
    plt.savefig('../Plots/an38/TimeSeries%02d%s_an38.pdf' % (ipar,parameter_shortname_list[ipar]),dpi=200)
    plt.close()



###########################################################################
###########################################################################
# I plot the integrated POC budget over time, for the layer 0-200m and 200-600m
###########################################################################
###########################################################################

POC_0_200=MiP_POC_0_200+MaP_POC_0_200+bbp_POC_0_200
POC_200_600=MiP_POC_200_600+MaP_POC_200_600+bbp_POC_200_600
POC_0_200_extended=MiP_POC_extended_0_200+MaP_POC_0_200+bbp_POC_0_200
POC_200_600_extended=MiP_POC_extended_200_600+MaP_POC_200_600+bbp_POC_200_600
POC_0_200[POC_0_200<0]=0
POC_200_600[POC_200_600<0]=0
POC_0_200_extended[POC_0_200_extended<0]=0
POC_200_600_extended[POC_200_600_extended<0]=0
# Parameters for the plot
width, height = 0.8, 0.5

#POC 0-200 vs 200-600
set_ylim_lower, set_ylim_upper = min(POC_0_200.min(),POC_200_600.min()),70 #max(POC_0_200.max(),POC_200_600.max())
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
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
ax.text(-0.075, 1.05, 'b', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an38/IntegratedPOC_vs_time_an38.pdf' ,dpi=200)
plt.close()

#POC 0-200 vs 200-600 extended
set_ylim_lower, set_ylim_upper = min(POC_0_200_extended.min(),POC_200_600_extended.min()),70 #max(POC_0_200_extended.max(),POC_200_600_extended.max())
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,POC_0_200_extended,'r',label='0-200 m')
plt.scatter(list_dates,POC_0_200_extended,c='r')
plt.plot(list_dates,POC_200_600_extended,'b',label='200-600 m\n(not magnified)')
plt.scatter(list_dates,POC_200_600_extended,c='b')
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
ax.text(-0.075, 1.05, 'b', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an38/IntegratedPOC_extended_vs_time_an38.pdf' ,dpi=200)
plt.close()



# 0-200 m layer
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,bbp_POC_0_200,'y',label='bbpPOC')
plt.scatter(list_dates,bbp_POC_0_200,c='y')
plt.plot(list_dates,MiP_POC_0_200,'c',label='MiP')
plt.scatter(list_dates,MiP_POC_0_200,c='c')
plt.plot(list_dates,MiP_POC_extended_0_200,'m',label='MiP extended')
plt.scatter(list_dates,MiP_POC_extended_0_200,c='m')
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
plt.savefig('../Plots/an38/0_200m_MiP_MaP_bbpPOC_vs_time_an38.pdf' ,dpi=200)
plt.close()




# 200-600 m layer
set_ylim_lower, set_ylim_upper = POC_200_600.min(),15# POC_200_600.max()
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
plt.plot(list_dates,bbp_POC_200_600,'y',label='bbpPOC')
plt.scatter(list_dates,bbp_POC_200_600,c='y')
plt.plot(list_dates,MiP_POC_200_600,'c',label='MiP')
plt.scatter(list_dates,MiP_POC_200_600,c='c')
plt.plot(list_dates,MiP_POC_extended_200_600,'m',label='MiP extended')
plt.scatter(list_dates,MiP_POC_extended_200_600,c='m')
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
plt.savefig('../Plots/an38/200_600m_MiP_MaP_bbpPOC_vs_time_an38.pdf' ,dpi=200)
plt.close()

