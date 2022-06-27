import numpy as np
import pandas as pd
import os,sys
import netCDF4 as nc
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename_coriolis='6903095_Sprof_all.nc'

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 01
########################################################################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 02
########################################################################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
######### Fig. 02a
########################################################################################################################
# filename1=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
# # I read the data
# data=pd.read_csv(filename1, sep=',', header=0)




########################################################################################################################
######### Fig. 02b,c,d,f
########################################################################################################################
# import datetime
# import seawater as sw
# import gsw
# from scipy.signal import savgol_filter
# from scipy.interpolate import griddata
#
# #Here I define the time at which I want to finish the time series in the plot
# day_end_timeseries=np.array([2021,9,24])
# day_end_timeseries=matlab_datenum(day_end_timeseries)
#
# #######################################################################
# # I load the Coriolis data
# #######################################################################
# ds = nc.Dataset('%s/%s' % (storedir,filename_coriolis))
# lon=np.array(ds.variables['LONGITUDE'])
# lat=np.array(ds.variables['LATITUDE'])
# Date_Num=np.array(ds.variables['JULD'])
# date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
#
# #Standard variables
# temp=np.array(ds.variables['TEMP_ADJUSTED'])
# pres=np.array(ds.variables['PRES_ADJUSTED'])
# psal=np.array(ds.variables['PSAL_ADJUSTED'])
#
# #BGC Variables
# chla=np.array(ds.variables['CHLA_ADJUSTED'])
# doxy=np.array(ds.variables['DOXY_ADJUSTED'])
# bbp700=np.array(ds.variables['BBP700_ADJUSTED'])
#
# #If adjusted values are not available yet, I take the non adjusted ones
# if np.sum(temp==99999)==temp.size:
#     print('Taking non adjusted temperature')
#     temp = np.array(ds.variables['TEMP'])
# if np.sum(pres==99999)==pres.size:
#     print('Taking non adjusted pressure')
#     pres = np.array(ds.variables['PRES'])
# if np.sum(psal==99999)==psal.size:
#     print('Taking non adjusted salinity')
#     psal = np.array(ds.variables['PSAL'])
# if np.sum(chla==99999)==chla.size:
#     print('Taking non adjusted chlorophyll-a')
#     chla = np.array(ds.variables['CHLA'])
# if np.sum(doxy==99999)==doxy.size:
#     print('Taking non adjusted oxygen')
#     doxy = np.array(ds.variables['DOXY'])
# if np.sum(bbp700==99999)==bbp700.size:
#     print('Taking non adjusted bbp700')
#     bbp700 = np.array(ds.variables['BBP700'])
#
# #######################################################################
# #I tranform the pressure to depth
# #######################################################################
# mask_depth=pres!=99999 #I select only valid values
# lat_tmp=np.tile(lat,[pres.shape[1],1]).T
# lat_tmp=lat_tmp[mask_depth]
# pres_tmp=pres[mask_depth]
# depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
# depth=np.ones(temp.shape)*99999
# depth[mask_depth]=depth_tmp
#
# #I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
# #salinity and temperature first
# mask_dens=np.logical_and(pres!=99999,temp!=99999,psal!=99999) # I exclude the points with value = 99999
# lat_tmp=np.tile(lat,[pres.shape[1],1]).T
# lon_tmp=np.tile(lon,[pres.shape[1],1]).T
# lat_tmp=lat_tmp[mask_dens]
# lon_tmp=lon_tmp[mask_dens]
# pres_tmp=pres[mask_dens]
# psal_tmp=psal[mask_dens]
# temp_tmp=temp[mask_dens]
# abs_psal_tmp = gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp)  # I compute absolute salinity
# cons_tmp = gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)          # I compute conservative temperature
# dens_tmp = gsw.density.sigma0(abs_psal_tmp, cons_tmp)
# abs_psal=np.ones(temp.shape)*99999
# abs_psal[mask_dens]=abs_psal_tmp
# cons_temp=np.ones(temp.shape)*99999
# cons_temp[mask_dens]=cons_tmp
# dens=np.ones(temp.shape)*99999
# dens[mask_dens]=dens_tmp+1000
#
# #######################################################################
# # I transform the bbp700 to small POC (sPOC)
# #######################################################################
# from oceanpy import bbp700toPOC
# sPOC=bbp700.copy()*0+99999
# i=0
# for i in range(0,bbp700.shape[0]):
#     bbp700tmp=bbp700[i,:]
#     depth_tmp=depth[i,:]
#     temp_tmp=temp[i,:]
#     # I exclude nan values
#     sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)
#     bbp700tmp=bbp700tmp[sel]
#     depth_tmp=depth_tmp[sel]
#     temp_tmp=temp_tmp[sel]
#     # I convert to small POC (sPOC) and I set to 0 values <0
#     sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
#     sPOC_tmp[sPOC_tmp<0]=0
#     sPOC[i,sel]=sPOC_tmp
#
# #######################################################################
# # I select the data only when the BGC Argo float was inside the eddy AND before day_end_timeseries (which fixes the x limit)
# #######################################################################
# filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
# data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)
#
# sel_insideEddy = data_dist_radius['sel_insideEddy']
# datenum_profiles = data_dist_radius['Datenum']
# sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)
#
# lon=lon[sel_insideEddy]
# lat=lat[sel_insideEddy]
# Date_Num=Date_Num[sel_insideEddy]
# pres=pres[sel_insideEddy]
# depth=depth[sel_insideEddy,:]
# dens=dens[sel_insideEddy,:]
# cons_temp=cons_temp[sel_insideEddy,:]
# chla=chla[sel_insideEddy,:]
# doxy=doxy[sel_insideEddy,:]
# sPOC=sPOC[sel_insideEddy,:]
#
# #######################################################################
# # I calculate the mixed layer depth
# #######################################################################
# from oceanpy import mixed_layer_depth
# mld=np.array([])
# i=0
# for i in range(0,chla.shape[0]):
#     depth_tmp=depth[i,:]
#     temp_tmp=cons_temp[i,:]
#     # I exclude nan values
#     sel_non_nan=(depth_tmp!=99999)&(temp_tmp!=99999)
#     temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
#     mld_tmp,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
#     mld=np.append(mld,mld_tmp)
#
# #######################################################################
# # I load the critical depth
# #######################################################################
# a_file = open("%s/an45/data_an45.pkl" % storedir, "rb")
# data_an45 = pickle.load(a_file)
# critical_depth=data_an45['critical_depth']
# critical_depth_datenum=data_an45['critical_depth_datenum']
# critical_depth_datenum = critical_depth_datenum[~np.isnan(critical_depth)]
# critical_depth = critical_depth[~np.isnan(critical_depth)]
#
# #######################################################################
# # I plot
# #######################################################################
# day_start_eddy_merging=np.array([2021,8,1])
# day_end_eddy_merging=np.array([2021,8,11])
# day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
# day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)
#
# parameter_ylabel_list=['Temperature ($^{\circ}$C)','Chlorophyll-a (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)']
# parameter_panellabel_list=['b','d','c','f']
# parameter_shortname_list=['cons_temp','chla','doxy','bbpPOC']
# ipar=0
# for ipar in range(0,parameter_ylabel_list.__len__()):
#     if ipar==0:   parameter=cons_temp.copy()
#     elif ipar == 1: parameter=chla.copy()
#     elif ipar == 2: parameter=doxy.copy()
#     elif ipar == 3: parameter=sPOC.copy()
#
#     #I filter the profiles
#     parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([])
#     i=0
#     for i in range(0,parameter.shape[0]):
#         z = parameter[i, :]
#         sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
#         if ipar == 3: sel = (sel) & (z <= 100)
#         if sum(sel) > 0:
#             z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
#             z = savgol_filter(z, 5, 1)
#             parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
#             depth_parameter = np.concatenate((depth_parameter, y1))
#
#     parameter_filtered[parameter_filtered<0]=0
#     # I define the x and y arrays for the contourf plot
#     x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
#     y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
#     # I interpolate
#     x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
#     parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
#
#
#     ########################################################
#     ####### I plot: versus depth
#     ########################################################
#     if ipar==3:
#         parameter_interp_depth[parameter_interp_depth > 40] = 40
#
#     width, height = 0.8, 0.7
#     set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
#     fig = plt.figure(1, figsize=(12,8))
#     ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
#     ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
#     if (ipar==0):
#         plt.plot(Date_Num,mld,'w')
#     elif ipar==1:
#         plt.plot(critical_depth_datenum,critical_depth,'w');plt.plot(critical_depth_datenum,critical_depth,'w.')
#
#     plt.gca().invert_yaxis()
#     plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
#     plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
#     # draw colorbar
#     cbar = plt.colorbar(plot2)
#     cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
#     plt.ylabel('Depth (m)', fontsize=18)
#     #plt.title('%smm' % NP_sizeclass, fontsize=18)
#     #I set xticks
#     nxticks=10
#     xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
#     xticklabels=[]
#     for i in xticks:
#         date_time_obj = date_reference + datetime.timedelta(days=i)
#         xticklabels.append(date_time_obj.strftime('%d %B'))
#     ax.set_xticks(xticks)
#     ax.set_xticklabels(xticklabels)
#     plt.xticks(rotation=90,fontsize=12)
#     # I add the panel label
#     ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
#     # I add the grid
#     plt.grid(color='k', linestyle='dashed', linewidth=0.5)
#     # I save
#     plt.savefig('../Plots/Fig_Main_v01/Fig02%s_v01.pdf' % parameter_panellabel_list[ipar],dpi=200)
#     plt.close()
#
#
#


########################################################################################################################
######### Fig. 02e,g,h
########################################################################################################################
import datetime
import seawater as sw
import gsw
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import calendar

#Here I define the time at which I want to finish the time series in the plot
# day_end_timeseries=datetime.datetime(2021,9,24)
# day_end_timeseries = calendar.timegm(day_end_timeseries.timetuple())
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)

#######################################################################
# I load the Coriolis data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename_coriolis))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num_bbp=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")

Date_Vec=np.zeros([Date_Num_bbp.size,6])
for i in range(0,Date_Num_bbp.size):
    date_time_obj = date_reference + datetime.timedelta(days=Date_Num_bbp[i]-matlab_datenum(1950,1,1))
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
#######################################################################
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth_bbp=np.ones(temp.shape)*99999
depth_bbp[mask_depth]=depth_tmp

#######################################################################
# I transform the bbp700 to small POC (sPOC)
#######################################################################
from oceanpy import bbp700toPOC
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
# I convert the bbp dates to float values (in seconds from 1970 1 1)
#######################################################################
Date_Num_bbp_calendar = Date_Num_bbp.copy()
for i in range(0, Date_Num_bbp_calendar.size):
    date_time_obj = datetime.datetime(Date_Vec[i, 0], Date_Vec[i, 1], Date_Vec[i, 2],
                             Date_Vec[i, 3], Date_Vec[i, 4], Date_Vec[i, 5])
    Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.utcfromtimestamp(Date_Num[i])

#######################################################################
# I load the MiP MaP data
#######################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data_ecopart=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data_ecopart.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data_ecopart
lon=np.array(data_ecopart['Longitude'][sel_filename])
lat=np.array(data_ecopart['Latitude'][sel_filename])
Date_Time=np.array(data_ecopart['Date_Time'][sel_filename])
pressure=np.array(data_ecopart['Pressure [dbar]'][sel_filename])
Flux=np.array(data_ecopart['Flux_mgC_m2'][sel_filename])
MiP_abund=np.array(data_ecopart['MiP_abun'][sel_filename])
MaP_abund=np.array(data_ecopart['MaP_abun'][sel_filename])
MiP_POC=np.array(data_ecopart['Mip_POC_cont_mgC_m3'][sel_filename])
MaP_POC=np.array(data_ecopart['Map_POC_cont_mgC_m3'][sel_filename])
depth=np.array(data_ecopart['Depth [m]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[sel_insideEddy[0:list_dates.size]]
Date_Num_bbp=Date_Num_bbp[sel_insideEddy]
Date_Num_bbp_calendar=Date_Num_bbp_calendar[sel_insideEddy]
depth_bbp=depth_bbp[sel_insideEddy]
temp=temp[sel_insideEddy]
bbp_POC=bbp_POC[sel_insideEddy,:]

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


POC_0_200_int=MiP_POC_0_200_int+MaP_POC_0_200_int+bbp_POC_0_200_int
POC_200_600_int=MiP_POC_200_600_int+MaP_POC_200_600_int+bbp_POC_200_600_int

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(POC_0_200_int.min(),POC_200_600_int.min()*10),max(POC_0_200_int.max(),POC_200_600_int.max()*10)

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
plt.savefig('../Plots/an66/IntegratedPOC_vs_time_an66.pdf' ,dpi=200)
plt.close()




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################




########################################################################################################################
######### Fig. 02e
########################################################################################################################







