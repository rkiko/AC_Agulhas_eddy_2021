
########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 01
########################################################################################################################
########################################################################################################################
########################################################################################################################
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
# region Fig. 02b,c,d,f
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
######
import datetime
import seawater as sw
import gsw
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to finish the time series in the plot
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)

#######################################################################
# I load the Coriolis data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename_coriolis))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")

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
# I select the data only when the BGC Argo float was inside the eddy AND before day_end_timeseries (which fixes the x limit)
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

lon=lon[sel_insideEddy]
lat=lat[sel_insideEddy]
Date_Num=Date_Num[sel_insideEddy]
pres=pres[sel_insideEddy]
depth=depth[sel_insideEddy,:]
dens=dens[sel_insideEddy,:]
cons_temp=cons_temp[sel_insideEddy,:]
chla=chla[sel_insideEddy,:]
doxy=doxy[sel_insideEddy,:]
sPOC=sPOC[sel_insideEddy,:]

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from oceanpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,chla.shape[0]):
    depth_tmp=depth[i,:]
    temp_tmp=cons_temp[i,:]
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
critical_depth_datenum=data_an45['critical_depth_datenum']
critical_depth_datenum = critical_depth_datenum[~np.isnan(critical_depth)]
critical_depth = critical_depth[~np.isnan(critical_depth)]

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)

parameter_ylabel_list=['Temperature ($^{\circ}$C)','Chlorophyll-a (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)']
parameter_panellabel_list=['b','d','c','f']
parameter_shortname_list=['cons_temp','chla','doxy','bbpPOC']
ipar=0
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0:   parameter=cons_temp.copy()
    elif ipar == 1: parameter=chla.copy()
    elif ipar == 2: parameter=doxy.copy()
    elif ipar == 3: parameter=sPOC.copy()

    #I filter the profiles
    parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        z = parameter[i, :]
        sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
        if ipar == 3: sel = (sel) & (z <= 100)
        if sum(sel) > 0:
            z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
            z = savgol_filter(z, 5, 1)
            parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
            depth_parameter = np.concatenate((depth_parameter, y1))

    parameter_filtered[parameter_filtered<0]=0
    # I define the x and y arrays for the contourf plot
    x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
    y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),50)
    # I interpolate
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
    parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")


    ########################################################
    ####### I plot: versus depth
    ########################################################
    if ipar==3:
        parameter_interp_depth[parameter_interp_depth > 40] = 40

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth)
    if (ipar==0):
        plt.plot(Date_Num,mld,'w')
    elif ipar==1:
        plt.plot(critical_depth_datenum,critical_depth,'w');plt.plot(critical_depth_datenum,critical_depth,'w.')

    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
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
    plt.savefig('../Plots/Fig_Main_v02/Fig02%s_v02.pdf' % parameter_panellabel_list[ipar],dpi=200)
    plt.close()
# endregion





########################################################################################################################
######### Fig. 02e,g,h
########################################################################################################################
# region Fig. 02e,g,h
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
###########
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

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

ipar=0
parameter_shortname_list=['MiP_POC','MaP_POC','bbpPOC']
parameter_panellabel_list=['g','h','f']
parameter_ylabel_list=['MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)','$b_{bp}$POC (mgC $m^{-3}$)']
max_parameter_list=np.array([2.15,0.30,40])
MiP_POC_0_200=np.array([]);MiP_POC_200_600=np.array([])
MaP_POC_0_200=np.array([]);MaP_POC_200_600=np.array([])
bbp_POC_0_200=np.array([]);bbp_POC_200_600=np.array([])
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar == 0: parameter=MiP_POC.copy()
    elif ipar == 1: parameter=MaP_POC.copy()
    elif ipar == 2: parameter=bbp_POC.copy()

    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
    if ipar == 2:
        i=0
        for i in range(0, bbp_POC.shape[0]):
            z=parameter[i,:];y=depth_bbp[i,:];x = Date_Num_bbp_calendar[i]
            z[z>100] = 99999
            sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
            sel3 = z == 0
            if sum(sel2) > 0:
                z = savgol_filter(z, 5, 1)
                z[sel3] = 0
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, np.tile(x,sum(sel2)) ))
                depth_filtered = np.concatenate((depth_filtered, y2))
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
                z=savgol_filter(z,5,1)
                parameter_filtered = np.concatenate((parameter_filtered, z))
                Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
                depth_filtered = np.concatenate((depth_filtered, y2))
                # sel_200 and sel_200_600 are used only for the POC integrated in time
                sel_0_200 = np.abs(y2) < 200
                sel_200_600 = (np.abs(y2) >= 200) & (np.abs(y2) <600)
                if ipar==0: MiP_POC_0_200=np.append(MiP_POC_0_200,np.mean(z[sel_0_200]));MiP_POC_200_600=np.append(MiP_POC_200_600,np.mean(z[sel_200_600]))
                if ipar==1: MaP_POC_0_200=np.append(MaP_POC_0_200,np.mean(z[sel_0_200]));MaP_POC_200_600=np.append(MaP_POC_200_600,np.mean(z[sel_200_600]))

    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    sel_0_200 = (np.abs(y_filtered) >= 0) & (np.abs(y_filtered) < 200)
    sel_200_600 = (np.abs(y_filtered) >= 200) & (np.abs(y_filtered) < 600)
    if ipar==0: MiP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MiP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==1: MaP_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);MaP_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)
    if ipar==2: bbp_POC_0_200_int = np.mean(parameter_interp[sel_0_200, :], 0);bbp_POC_200_600_int = np.mean(parameter_interp[sel_200_600, :], 0)

    if ipar == 2: continue

    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
    plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
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
        xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90,fontsize=12)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/Fig_Main_v02/Fig02%s_v02.pdf' % (parameter_panellabel_list[ipar]),dpi=200)
    plt.close()


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
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
ax.text(-0.075, 1.05, 'e', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v02/Fig02e_v02.pdf' ,dpi=200)
plt.close()
# endregion




########################################################################################################################
######### Fig. 02e_v01b
########################################################################################################################
# region Fig. 02e_v01b

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
###########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to start and end the time series in the plot
day_start_timeseries=np.array([2021,4,13])
day_start_timeseries=matlab_datenum(day_start_timeseries)
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
ndays=int(day_end_timeseries-day_start_timeseries+1)

#######################################################################
# I load the mixed layer depth
#######################################################################
dictionary_data = {"mld": mld,"Date_Num": Date_Num,"lon": lon,"lat": lat}
a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an68/data_MLD_an68.pkl" % (home), "rb")

#######################################################################
# I load and process MiP and MaP data
#######################################################################

#I load the file with the flux and POC
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
Date_Time=np.array(data['Date_Time'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])
dens=np.array(data['Potential density [kg/m3]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
filename_dist_radius=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()
data_dist_radius=pd.read_csv(filename_dist_radius, sep=',', header=0)

sel_insideEddy = data_dist_radius['sel_insideEddy']
datenum_profiles = data_dist_radius['Datenum']
sel_insideEddy = (datenum_profiles<=day_end_timeseries)&(sel_insideEddy==1)

list_dates=list_dates[sel_insideEddy[0:list_dates.size]]

n_profiles=list_dates.size

#######################################################################
# I load and process bbp data
#######################################################################

# I load the bbp data and I select only those in the prescribed period
storedir = '%s/GIT/AC_Agulhas_eddy_2021/Data' % home
a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
data_an18 = pickle.load(a_file)
bbp_POC = data_an18['bbp_POC']
Date_Num_bbp = data_an18['Date_Num_bbp']
Date_Vec_bbp = data_an18['Date_Vec_bbp']
depth_bbp = data_an18['depth_bbp']
dens_bbp = data_an18['dens_bbp']
a_file.close()

sel_dates = sel_insideEddy[0:Date_Num_bbp.size]
Date_Num_bbp = Date_Num_bbp[sel_dates]
Date_Vec_bbp = Date_Vec_bbp[sel_dates,:]
depth_bbp = depth_bbp[sel_dates,:]
dens_bbp = dens_bbp[sel_dates,:]
bbp_POC = bbp_POC[sel_dates, :]

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num_bbp_calendar = Date_Num_bbp.copy()
for i in range(0,Date_Num_bbp_calendar.size):
    date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                             Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
    Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.utcfromtimestamp(Date_Num[i])

########################################################################################################################
# Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
# interpolate it over a regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because
# some profiles only reach 400 m; (iii) for each density value of the density grid, I calculate the corresponding depth.
# Finally, (iv) I extract the mean MiP+MaP+bbp values between depth0 and depthf and between day0 and dayf (I obtain a
# time series)
########################################################################################################################

##############################################
# Step 1 and 2, filter and interpolation
MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])

i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=MiP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        MiP_filtered = np.concatenate((MiP_filtered, z))
        Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x2))
        dens_MiP_filtered = np.concatenate((dens_MiP_filtered, y2))
    z=MaP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        MaP_filtered = np.concatenate((MaP_filtered, z))
        Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x2))
        dens_MaP_filtered = np.concatenate((dens_MaP_filtered, y2))

i=0
for i in range(0, bbp_POC.shape[0]):
    z=bbp_POC[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
    z[z>100] = 99999
    sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
    sel3=z==0
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        z[sel3]=0
        bbp_filtered = np.concatenate((bbp_filtered, z))
        Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, np.tile(x,sum(sel2)) ))
        dens_bbp_filtered = np.concatenate((dens_bbp_filtered, y2))

# I define the x and y arrays for the MiP+MaP+bbp interpolation
x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), 100)
y_filtered = np.linspace(dens_bbp_filtered.min(), dens_MaP_filtered.max(), 200)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
MiP_interp = griddata((Date_Num_MiP_filtered, dens_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
MaP_interp = griddata((Date_Num_MaP_filtered, dens_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
bbp_interp = griddata((Date_Num_bbp_filtered, dens_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")

##############################################
# Step 3: for each density layer I calculate the corresponding mean depth
sel = (Date_Num <= list_dates.max())&(Date_Num > list_dates.min())
dens_tmp=dens[sel];depth_tmp=depth[sel]
sel2=(~np.isnan(dens_tmp))&(~np.isnan(depth_tmp))&(dens_tmp!=99999)&(depth_tmp!=99999)
dens_tmp=dens_tmp[sel2];depth_tmp=depth_tmp[sel2]

sel_bbp=(~np.isnan(dens_bbp))&(~np.isnan(depth_bbp))&(dens_bbp!=99999)&(depth_bbp!=99999)
dens_tmp_bbp=dens_bbp[sel_bbp];depth_tmp_bbp=depth_bbp[sel_bbp]

depth_bbp_filtered=np.array([]);depth_MiP_filtered=np.array([])
i=0
for i in range(0,y_filtered.size):
    if i==0:    dens0=y_filtered[i]
    else:       dens0=abs( y_filtered[i]+y_filtered[i-1] )*0.5
    if i==(y_filtered.size-1):  dens1=y_filtered[i]
    else:                       dens1=abs( y_filtered[i]+y_filtered[i+1] )*0.5
    #Depth calculation for MiP and MaP
    sel_dens = (dens_tmp >= dens0) & (dens_tmp < dens1)
    if sum(sel_dens) > 0:
        depth_MiP_filtered = np.append(depth_MiP_filtered, np.nanmean(depth_tmp[sel_dens]))
    else:
        depth_MiP_filtered = np.append(depth_MiP_filtered, np.array([np.nan]))

    # Depth calculation for bbp
    sel_dens = (dens_tmp_bbp >= dens0) & (dens_tmp_bbp < dens1)
    if sum(sel_dens) > 0:
        depth_bbp_filtered = np.append(depth_bbp_filtered, np.nanmean(depth_tmp_bbp[sel_dens]))
    else:
        depth_bbp_filtered = np.append(depth_bbp_filtered, np.array([np.nan]))



##############################################
# Step 4, I calculate the mean MiP+MaP+bbp (and std) between depth0 and depthf between day0 and dayf
sel_0_200 = (np.abs(depth_MiP_filtered) >= 0) & (np.abs(depth_MiP_filtered) < 200)
sel_0_200_bbp = (np.abs(depth_bbp_filtered) >= 0) & (np.abs(depth_bbp_filtered) < 200)
MiP_POC_0_200=np.mean(MiP_interp[sel_0_200,:],0)
MaP_POC_0_200=np.mean(MaP_interp[sel_0_200,:],0)
bbp_POC_0_200=np.mean(bbp_interp[sel_0_200_bbp,:],0)

MiP_POC_0_200_std = np.std(MiP_interp[sel_0_200, :], 0)
MaP_POC_0_200_std = np.std(MaP_interp[sel_0_200, :], 0)
bbp_POC_0_200_std = np.std(bbp_interp[sel_0_200_bbp, :], 0)

sel_200_600 = (np.abs(depth_MiP_filtered) >= 200) & (np.abs(depth_MiP_filtered) < 600)
sel_200_600_bbp = (np.abs(depth_bbp_filtered) >= 200) & (np.abs(depth_bbp_filtered) < 600)
MiP_POC_200_600=np.mean(MiP_interp[sel_200_600,:],0)
MaP_POC_200_600=np.mean(MaP_interp[sel_200_600,:],0)
bbp_POC_200_600=np.mean(bbp_interp[sel_200_600_bbp,:],0)

MiP_POC_200_600_std = np.std(MiP_interp[sel_200_600, :], 0)
MaP_POC_200_600_std = np.std(MaP_interp[sel_200_600, :], 0)
bbp_POC_200_600_std = np.std(bbp_interp[sel_200_600_bbp, :], 0)

Integrated_POC_mgC_m3_0_200 = MiP_POC_0_200 + MaP_POC_0_200 + bbp_POC_0_200
Integrated_POC_mgC_m3_0_200_std = np.sqrt( MiP_POC_0_200_std**2 + MaP_POC_0_200_std**2 + bbp_POC_0_200_std**2 )
Integrated_POC_mgC_m3_200_600 = MiP_POC_200_600 + MaP_POC_200_600 + bbp_POC_200_600
Integrated_POC_mgC_m3_200_600_std = np.sqrt( MiP_POC_200_600_std**2 + MaP_POC_200_600_std**2 + bbp_POC_200_600_std**2 )
list_dates_Integrated_POC = x_filtered.copy()

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = min(Integrated_POC_mgC_m3_0_200.min(),Integrated_POC_mgC_m3_200_600.min()*10),max(Integrated_POC_mgC_m3_0_200.max(),Integrated_POC_mgC_m3_200_600.max()*10)

fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, set_ylim_upper*1.1), xlim=(list_dates.min(), list_dates.max()))
# plt.plot(list_dates,POC_0_200,'r',label='0-200 m')
# plt.scatter(list_dates,POC_0_200,c='r')
# plt.plot(list_dates,POC_200_600*10,'b',label='200-600 m')
# plt.scatter(list_dates,POC_200_600*10,c='b')
plt.plot(list_dates_Integrated_POC,Integrated_POC_mgC_m3_0_200,'r',linewidth=3,label='0-200 m')
plt.plot(list_dates_Integrated_POC,Integrated_POC_mgC_m3_200_600*10,'b',linewidth=3,label='200-600 m')
plt.vlines(day_start_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
plt.vlines(day_end_eddy_merging, ymin=0, ymax=600, color='k',linestyles='dashed',linewidth=3)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=14)
plt.ylabel('Average POC (mgC/m$^3$)', fontsize=15)
ax.text(-0.075, 1.05, 'e', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v02/Fig02e_v01B.pdf' ,dpi=200)
plt.close()



# endregion

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 03
########################################################################################################################
########################################################################################################################
########################################################################################################################

#######################################################################################################################
######### Fig. 03a,b
########################################################################################################################
# region Fig. 03a,b
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
#######
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

#Here I define the time at which I want to finish the time series in the plot
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
delta_depth=15                  # around of the depth which I consider when extracting the flux

#######################################################################
# I load the ecopart data
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

# I extract the data
lon=np.array(data_ecopart['Longitude'][sel_filename])
lat=np.array(data_ecopart['Latitude'][sel_filename])
Date_Time=np.array(data_ecopart['Date_Time'][sel_filename])
Flux=np.array(data_ecopart['Flux_mgC_m2'][sel_filename])
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

#######################################################################
# I filter and interpolate the flux
#######################################################################

day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

parameter=Flux.copy()
parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
# I filter the flux prophiles
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
    z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        parameter_filtered = np.concatenate((parameter_filtered, z))
        Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
        depth_filtered = np.concatenate((depth_filtered, y2))

x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
# I interpolate
parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

depthf = 600  # final depth
sel_depthf_600 = (np.abs(y_filtered) > depthf - delta_depth) & (np.abs(y_filtered) <= depthf + delta_depth)
Flux_interp_600 = np.mean(parameter_interp[sel_depthf_600, :], axis=0)

depthf = 200  # final depth
sel_depthf_200 = (np.abs(y_filtered) > depthf - delta_depth) & (np.abs(y_filtered) <= depthf + delta_depth)
Flux_interp_200 = np.mean(parameter_interp[sel_depthf_200, :], axis=0)

########################################################################################################################
######### Fig. 03a
########################################################################################################################
fs=9
width, height = 0.82, 0.8

# First plot: flux calculated without considering smallest size classes and with old eta and b values
fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x_filtered,Flux_interp_200,'r',label='Flux at 200 m')
plt.plot(x_filtered,Flux_interp_600,'b',label='Flux at 600 m')
plt.xlim(x_filtered.min(),x_filtered.max())
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
ax.text(-0.115, 1.05, 'a', transform=ax.transAxes,fontsize=14, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.ylabel('Flux(mgC/$m^2$/d)',fontsize=7)
plt.legend(fontsize=7)#,ncol=2)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
ax.set_xticks(xticks)
ax.set_xticklabels([])
plt.xticks(rotation=90, fontsize=7)
plt.savefig('../Plots/Fig_Main_v02/Fig03a_v02.pdf'  ,dpi=200)
plt.close()


########################################################################################################################
######### Fig. 03b
########################################################################################################################
max_Flux=32

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
parameter_plot=parameter_interp.copy()
parameter_plot[parameter_plot<0]=0
parameter_plot[parameter_plot>max_Flux]=max_Flux
ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
plt.gca().invert_yaxis()
plt.vlines(day_start_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
plt.vlines(day_end_eddy_merging,ymin=0,ymax=600,color='w',linestyles='dashed')
# I draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the panel label
ax.text(-0.05, 1.05, 'b', transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v02/Fig03b_v02.pdf' ,dpi=200)
plt.close()
# endregion

#######################################################################################################################
######### Fig. 03a_v01B
########################################################################################################################
#region Fig. 03a_v01B
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
###########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
delta_depth_flux = 15  # around of the depth which I consider when extracting the flux

#Here I define the time at which I want to start and end the time series in the plot
day_start_timeseries=np.array([2021,4,13])
day_start_timeseries=matlab_datenum(day_start_timeseries)
day_end_timeseries=np.array([2021,9,24])
day_end_timeseries=matlab_datenum(day_end_timeseries)
ndays=int(day_end_timeseries-day_start_timeseries+1)

#######################################################################
# I load and process MiP and MaP data
#######################################################################

#I load the file with the flux and POC
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
Date_Time=np.array(data['Date_Time'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])
dens=np.array(data['Potential density [kg/m3]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

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

n_profiles=list_dates.size

########################################################################################################################
# Here I extract the flux values at depth0 and depthf. To do so, (i) I filter it with a savgol function, then (ii) I
# interpolate it over a regular grid in time and density. This step is necessary to have the flux at 600 m, because some
# profiles only reach 400 m; (iii) for each density value of the density grid, I calculate the corresponding depth.
# Finally, (iv) I extract the flux values at depth0 and depthf
########################################################################################################################

##############################################
# Step 1 and 2, filter and interpolation
Flux_filtered=np.array([]);dens_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=dens[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        Flux_filtered = np.concatenate((Flux_filtered, z))
        Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
        dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y2))

# I define the x and y arrays for the Flux interpolation
x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
y_filtered = np.linspace(dens_Flux_filtered.min(), dens_Flux_filtered.max(), 200)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
Flux_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")

##############################################
# Step 3: for each density layer I calculate the corresponding mean depth
depth_Flux_filtered=np.array([])
i=0
for i in range(0,y_filtered.size):
    if i==0:    dens0=y_filtered[i]
    else:       dens0=abs( y_filtered[i]+y_filtered[i-1] )*0.5
    if i==(y_filtered.size-1):  dens1=y_filtered[i]
    else:                       dens1=abs( y_filtered[i]+y_filtered[i+1] )*0.5
    #Depth calculation
    depth_Flux_filtered_tmp = np.array([])
    j=0
    for j in range(0,list_dates.size):
        sel=Date_Num==list_dates[j]
        dens_tmp=dens[sel];depth_tmp=depth[sel];sel2=(~np.isnan(dens_tmp))&(~np.isnan(depth_tmp));dens_tmp=dens_tmp[sel2];depth_tmp=depth_tmp[sel2]
        if sum(sel2) > 0:
            sel_dens = (dens_tmp >= dens0) & (dens_tmp < dens1)
            if sum(sel_dens) > 0:
                depth_Flux_filtered_tmp = np.append(depth_Flux_filtered_tmp, np.mean(depth_tmp[sel_dens]) )
            else:
                depth_Flux_filtered_tmp = np.append( depth_Flux_filtered_tmp, np.array([np.nan]) )
        else:
            depth_Flux_filtered_tmp = np.append( depth_Flux_filtered_tmp, np.array([np.nan]) )

    if sum(~np.isnan(depth_Flux_filtered_tmp))==0:
        depth_Flux_filtered = np.append(depth_Flux_filtered, np.array([np.nan]) )
    else:
        depth_Flux_filtered = np.append(depth_Flux_filtered, np.nanmean(depth_Flux_filtered_tmp))

##############################################
# Step 4, flux extraction at depth0 and depthf
depth0=200
depthf=600
sel_layer = (np.abs(depth_Flux_filtered) >= depth0-delta_depth_flux) & (np.abs(depth_Flux_filtered) < depth0+delta_depth_flux)
Flux_depth0 = np.mean(Flux_interp[sel_layer,:],axis=0)
sel_layer = (np.abs(depth_Flux_filtered) >= depthf - delta_depth_flux) & (np.abs(depth_Flux_filtered) < depthf + delta_depth_flux)
Flux_depthf = np.mean(Flux_interp[sel_layer,:],axis=0)


#######################################################################
# I plot
#######################################################################
day_start_eddy_merging = datetime.datetime(2021,8,1)
day_start_eddy_merging = calendar.timegm(day_start_eddy_merging.timetuple())
day_end_eddy_merging = datetime.datetime(2021,8,11)
day_end_eddy_merging = calendar.timegm(day_end_eddy_merging.timetuple())

fs=9
width, height = 0.82, 0.8

fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x_filtered,Flux_depth0,'r',label='Flux at 200 m')
plt.plot(x_filtered,Flux_depthf,'b',label='Flux at 600 m')
plt.xlim(x_filtered.min(),x_filtered.max())
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k',linestyles='dashed')
ax.text(-0.115, 1.05, 'a', transform=ax.transAxes,fontsize=14, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.ylabel('Flux(mgC/$m^2$/d)',fontsize=7)
plt.legend(fontsize=7)#,ncol=2)
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates.min(), list_dates.max(), nxticks)
xticklabels = []
ax.set_xticks(xticks)
ax.set_xticklabels([])
plt.xticks(rotation=90, fontsize=7)
plt.savefig('../Plots/Fig_Main_v02/Fig03a_v01B.pdf'  ,dpi=200)
plt.close()

#endregion


########################################################################################################################
########################################################################################################################
########################################################################################################################
######### FIG 04
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Fig. 04
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
########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import seawater as sw
import gsw
from lin_fit import lin_fit

#######################################################################
# I define the function for the carbon budget calculation
#######################################################################
def carbon_budget_calculation(dens0,densf,day0,dayf):
    ########################################################################################################################
    # Starting parameters
    ########################################################################################################################
    # dayf = day0+timedelta(days=ndays) # final date for the carbon budget calculation
    ndays = (dayf - day0).days  # number of days
    layer_thickness = densf - dens0
    delta_dens_flux = 0.015     # around of the density which I consider when extracting the flux
    Oxy2C = 0.89                # to convert from mol of oxygen to mol of carbon
    mol2gC = 12.0107            # to convert from mol of carbon to grams of carbon
    day0_float = calendar.timegm(day0.timetuple())
    dayf_float = calendar.timegm(dayf.timetuple())
    day0_datenum = matlab_datenum(day0.year,day0.month,day0.day,day0.hour,day0.minute,day0.second)
    dayf_datenum = matlab_datenum(dayf.year,dayf.month,dayf.day,dayf.hour,dayf.minute,dayf.second)
    delta_rho=0.025  # around of the isopycnal
    reference_isopycnal_list=np.r_[1026.3:1027.50001:delta_rho]

    ########################################################################################################################
    # I load and process data
    ########################################################################################################################

    #I load the file with the flux and POC
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
    Date_Time=np.array(data['Date_Time'][sel_filename])
    depth=np.array(data['Depth [m]'][sel_filename])
    dens=np.array(data['Potential density [kg/m3]'][sel_filename])
    Flux=np.array(data['Flux_mgC_m2'][sel_filename])
    Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
    Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
    Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
    MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
    MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
    MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num=np.r_[0:Flux.size]
    for i in Date_Num:
        date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    # I select the data only in the prescribed period
    list_dates=np.sort(np.unique(Date_Num))
    list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
    n_profiles=list_dates.size

    # I load the bbp data and I select only those in the prescribed period
    storedir = '%s/GIT/AC_Agulhas_eddy_2021/Data' % home
    a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
    data_an18 = pickle.load(a_file)
    bbp_POC = data_an18['bbp_POC']
    Date_Num_bbp = data_an18['Date_Num_bbp']
    Date_Vec_bbp = data_an18['Date_Vec_bbp']
    depth_bbp = data_an18['depth_bbp']
    dens_bbp = data_an18['dens_bbp']
    a_file.close()

    sel_dates = (Date_Num_bbp>=day0_datenum)&(Date_Num_bbp<=dayf_datenum)
    Date_Num_bbp = Date_Num_bbp[sel_dates]
    Date_Vec_bbp = Date_Vec_bbp[sel_dates,:]
    depth_bbp = depth_bbp[sel_dates,:]
    dens_bbp = dens_bbp[sel_dates,:]
    bbp_POC = bbp_POC[sel_dates, :]

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_bbp_calendar = Date_Num_bbp.copy()
    for i in range(0,Date_Num_bbp_calendar.size):
        date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                                 Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
        Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
        # datetime.utcfromtimestamp(Date_Num[i])

    ########################################################################################################################
    # Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because
    # some profiles only reach 400 m; (iii) I extract the mean MiP+MaP+bbp values between dens0 and densf and between day0 and
    # dayf (I obtain a time series)
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
    MiP_extended_filtered=np.array([]);dens_MiP_extended_filtered=np.array([]);Date_Num_MiP_extended_filtered=np.array([])
    MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
    bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])

    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=MiP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_filtered = np.concatenate((MiP_filtered, z))
            Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x2))
            dens_MiP_filtered = np.concatenate((dens_MiP_filtered, y2))
        z=MiP_POC_extended[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_extended_filtered = np.concatenate((MiP_extended_filtered, z))
            Date_Num_MiP_extended_filtered = np.concatenate((Date_Num_MiP_extended_filtered, x2))
            dens_MiP_extended_filtered = np.concatenate((dens_MiP_extended_filtered, y2))
        z=MaP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MaP_filtered = np.concatenate((MaP_filtered, z))
            Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x2))
            dens_MaP_filtered = np.concatenate((dens_MaP_filtered, y2))

    i=0
    for i in range(0, bbp_POC.shape[0]):
        z=bbp_POC[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>100] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_filtered = np.concatenate((bbp_filtered, z))
            Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_filtered = np.concatenate((dens_bbp_filtered, y2))

    # I define the x and y arrays for the MiP+MaP+bbp interpolation
    x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), ndays)
    y_filtered = np.linspace(dens_bbp_filtered.min(), dens_MaP_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    MiP_interp = griddata((Date_Num_MiP_filtered, dens_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MiP_extended_interp = griddata((Date_Num_MiP_extended_filtered, dens_MiP_extended_filtered), MiP_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MaP_interp = griddata((Date_Num_MaP_filtered, dens_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_interp = griddata((Date_Num_bbp_filtered, dens_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, I calculate the mean MiP+MaP+bbp (and std) between dens0 and densf between day0 and dayf
    sel_dens0_densf = (np.abs(y_filtered) >= dens0) & (np.abs(y_filtered) < densf)
    MiP_POC_dens0_densf=np.mean(MiP_interp[sel_dens0_densf,:],0)
    MiP_POC_extended_dens0_densf=np.mean(MiP_extended_interp[sel_dens0_densf,:],0)
    MaP_POC_dens0_densf=np.mean(MaP_interp[sel_dens0_densf,:],0)
    bbp_POC_dens0_densf=np.mean(bbp_interp[sel_dens0_densf,:],0)

    MiP_POC_dens0_densf_std = np.std(MiP_interp[sel_dens0_densf, :], 0)
    MiP_POC_extended_dens0_densf_std = np.std(MiP_extended_interp[sel_dens0_densf, :], 0)
    MaP_POC_dens0_densf_std = np.std(MaP_interp[sel_dens0_densf, :], 0)
    bbp_POC_dens0_densf_std = np.std(bbp_interp[sel_dens0_densf, :], 0)

    Integrated_POC_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    Integrated_POC_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    list_dates_Integrated_POC = x_filtered.copy()

    ########################################################################################################################
    # Here I extract the flux values at dens0 and densf. To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid in time and density. This step is necessary to have the flux at 600 m, because some
    # profiles only reach 400 m; (iii) I extract the flux values at dens0 and densf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    Flux_filtered=np.array([]);dens_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
    Flux_eta_b_filtered=np.array([]);dens_Flux_eta_b_filtered=np.array([]);Date_Num_Flux_eta_b_filtered=np.array([])
    Flux_extended_filtered=np.array([]);dens_Flux_extended_filtered=np.array([]);Date_Num_Flux_extended_filtered=np.array([])
    Flux_extended_filtered_eta_b=np.array([]);dens_Flux_extended_filtered_eta_b=np.array([]);Date_Num_Flux_extended_filtered_eta_b=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=Flux[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_filtered = np.concatenate((Flux_filtered, z))
            Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
            dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y2))
        z=Flux_eta_b[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_eta_b_filtered = np.concatenate((Flux_eta_b_filtered, z))
            Date_Num_Flux_eta_b_filtered = np.concatenate((Date_Num_Flux_eta_b_filtered, x2))
            dens_Flux_eta_b_filtered = np.concatenate((dens_Flux_eta_b_filtered, y2))
        z=Flux_extended[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered = np.concatenate((Flux_extended_filtered, z))
            Date_Num_Flux_extended_filtered = np.concatenate((Date_Num_Flux_extended_filtered, x2))
            dens_Flux_extended_filtered = np.concatenate((dens_Flux_extended_filtered, y2))
        z=Flux_extended_eta_b[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered_eta_b = np.concatenate((Flux_extended_filtered_eta_b, z))
            Date_Num_Flux_extended_filtered_eta_b = np.concatenate((Date_Num_Flux_extended_filtered_eta_b, x2))
            dens_Flux_extended_filtered_eta_b = np.concatenate((dens_Flux_extended_filtered_eta_b, y2))

    # I define the x and y arrays for the Flux interpolation
    x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
    y_filtered = np.linspace(dens_Flux_filtered.min(), dens_Flux_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    Flux_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_eta_b_interp = griddata((Date_Num_Flux_eta_b_filtered, dens_Flux_eta_b_filtered), Flux_eta_b_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_interp = griddata((Date_Num_Flux_extended_filtered, dens_Flux_extended_filtered), Flux_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_eta_b_interp = griddata((Date_Num_Flux_extended_filtered_eta_b, dens_Flux_extended_filtered_eta_b), Flux_extended_filtered_eta_b,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, flux extraction at dens0 and densf

    sel_layer = (np.abs(y_filtered) >= dens0-delta_dens_flux) & (np.abs(y_filtered) < dens0+delta_dens_flux)
    Flux_dens0 = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_eta_b_dens0 = np.mean(Flux_eta_b_interp[sel_layer, :], axis=0)
    Flux_extended_dens0 = np.mean(Flux_extended_interp[sel_layer, :], axis=0)
    Flux_extended_eta_b_dens0 = np.mean(Flux_extended_eta_b_interp[sel_layer, :], axis=0)

    sel_layer = (np.abs(y_filtered) >= densf - delta_dens_flux) & (np.abs(y_filtered) < densf + delta_dens_flux)
    Flux_densf = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_eta_b_densf = np.mean(Flux_eta_b_interp[sel_layer,:],axis=0)
    Flux_extended_densf = np.mean(Flux_extended_interp[sel_layer,:],axis=0)
    Flux_extended_eta_b_densf = np.mean(Flux_extended_eta_b_interp[sel_layer,:],axis=0)


    ########################################################################################################################
    # Here I calculate the carbon consumption rate due to (i) oxygen consumption and (ii) PARR
    ########################################################################################################################

    ############### I load Coriolis data with the oxygen information
    filename='6903095_Sprof_all.nc'
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

    ############### I load the PARR correlated data (density, depth, and date) from Ecopart
    Date_Time_PARR=np.array(data['Date_Time'][sel_filename])
    depth_PARR=np.array(data['Depth [m]'][sel_filename])
    dens_PARR=np.array(data['Potential density [kg/m3]'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_PARR=np.r_[0:Date_Time_PARR.size]
    for i in range(0,Date_Time_PARR.size):
        date_time_obj = datetime.datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    list_dates_PARR=np.unique(Date_Num_PARR)

    #############################################
    ############### Loop on the different isopycnal values chosen for the study of the oxygen profile
    Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+ndays]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))
    depth_isopycnal = np.array([])
    slopes_list_doxy = np.array([])
    slopes_ci_list_doxy = np.array([])
    reference_isopycnal_list_sel = np.array([])

    i_isop=0
    for i_isop in range(0,reference_isopycnal_list.size):
        #Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
        #isopycnal, given by reference_isopycnal, and then the PARR concentration at that isopycnal
        reference_isopycnal=reference_isopycnal_list[i_isop]
        reference_isopycnal_down=reference_isopycnal-delta_rho/2
        if reference_isopycnal_down<reference_isopycnal_list[0]:  reference_isopycnal_down=reference_isopycnal_list[0]
        reference_isopycnal_up=reference_isopycnal+delta_rho/2
        if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]
        doxy_isopycnal=np.array([]);depth_isopycnal_tmp=np.array([]);Date_Num_isopycnal=np.array([])
        i=0
        for i in range(0,doxy.shape[0]):
            #Here, for the i-th profile, I select the oxygen, density and depth profiles of Coriolis data, excluding the nan values
            sel = (doxy[i, :] != 99999) & (dens[i, :] != 99999)
            z=doxy[i,sel];y=dens[i,sel];d=depth[i,sel]

            # Here I proceed only if the date is inside the Date_Num_limit fixed
            if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

                # Here I extract the oxygen along the isopycnal
                sel_layer = (y >= reference_isopycnal_down) & (y < reference_isopycnal_up)
                if np.sum(sel_layer) > 0:  # If sel_layer has some True values, then I take as the doxy of this isopycnal the mean of the doxy values in correspondence with these layers
                    doxy_tmp = np.mean(z[sel_layer])
                    depth_isopycnal_tmp2 = np.mean(d[sel_layer])
                    doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                    Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                    depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, depth_isopycnal_tmp2)
                else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the oxygen concentration at that depth
                    depth_isopycnal_tmp2 = np.array([])
                    for iy in range(0, y.size - 1):
                        if y[iy] <= reference_isopycnal < y[iy + 1]:
                            dist = (reference_isopycnal - y[iy]) / (y[iy + 1] - y[iy])
                            doxy_tmp = z[iy] + (z[iy + 1] - z[iy]) * dist
                            d_tmp = d[iy] + (d[iy + 1] - d[iy]) * dist
                            doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                            depth_isopycnal_tmp2 = np.append(depth_isopycnal_tmp2, d_tmp)
                            Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                    if depth_isopycnal_tmp2.size > 0:
                        depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, np.mean(depth_isopycnal_tmp2))

        # If I have at least three points I interpolate them, so that the slope gives me the respiration rate
        # (in micromol/kg per day)
        if Date_Num_isopycnal.size>2:
            (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)
            # fig = plt.figure(1, figsize=(12, 8))
            # plot1 = plt.scatter(Date_Num_isopycnal,doxy_isopycnal)
            depth_isopycnal = np.append(depth_isopycnal ,np.mean(depth_isopycnal_tmp))
            slopes_list_doxy = np.append(slopes_list_doxy ,interpol.slope)
            reference_isopycnal_list_sel = np.append(reference_isopycnal_list_sel ,reference_isopycnal)
            if i_isop==0:
                slopes_ci_list_doxy = np.reshape(slpe_ci.copy(),(1,2))
            else:
                slopes_ci_list_doxy = np.concatenate( ( slopes_ci_list_doxy ,np.reshape(slpe_ci,(1,2)) ), axis=0)

    O2_resp_mgC_m2_vs_depth=slopes_list_doxy.copy()*reference_isopycnal_list_sel*Oxy2C*mol2gC*layer_thickness*ndays/1000
    O2_resp_mgC_m2_vs_depth_ci=slopes_ci_list_doxy.copy()*( np.tile(reference_isopycnal_list_sel,(2,1)).T) *Oxy2C*mol2gC*layer_thickness*ndays/1000

    #############################################
    ############### Loop on different respiration types used to estimate PARR
    #List of the different Respiration types present in data
    list_Respi_types = [match for match in data.columns if "Respi" in match]
    nRespi= len(list_Respi_types)  #number of respiration types

    POC_resp_mgC_m2_vs_depth_list=np.zeros((reference_isopycnal_list.size,nRespi))
    POC_resp_mgC_m2_std_vs_depth_list = np.zeros((reference_isopycnal_list.size, nRespi))

    iRespi=0
    for iRespi in range(0,nRespi):
        PARR_nmol_l_h = np.array(data[list_Respi_types[iRespi]][sel_filename])
        # I convert the PARR measured in micromol/kg/day
        PARR_micromol_kg_day = PARR_nmol_l_h.copy() / 1000 * 24 / (dens_PARR / 1000)

        #############################################
        ############### Loop on the different isopycnal values chosen for the study of the PARR profile
        list_depth_PARR = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        list_dens_PARR = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        PARR_list = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        PARR_std_list = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))

        i_isop=0
        for i_isop in range(0,reference_isopycnal_list.size):
            #Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
            #isopycnal, given by reference_isopycnal, and then the PARR concentration at that isopycnal
            reference_isopycnal=reference_isopycnal_list[i_isop]
            reference_isopycnal_down=reference_isopycnal-delta_rho/2
            if reference_isopycnal_down<reference_isopycnal_list[0]:  reference_isopycnal_down=reference_isopycnal_list[0]
            reference_isopycnal_up=reference_isopycnal+delta_rho/2
            if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]
            PARR_isopycnal=np.array([]);depth_PARR_tmp=np.array([]);dens_PARR_tmp=np.array([])
            i=0
            for i in range(0,list_dates_PARR.size):
                # Here, for the i-th profile, I select the PARR, density and depth profiles of Ecopart data, excluding the nan values. I do the same for the bbp
                sel_PARR=Date_Num_PARR==list_dates_PARR[i]
                z_PARR=PARR_micromol_kg_day[sel_PARR];y_PARR=dens_PARR[sel_PARR];d_PARR=depth_PARR[sel_PARR]
                sel_PARR = (~np.isnan(z_PARR)) & (~np.isnan(y_PARR))
                z_PARR = z_PARR[sel_PARR];y_PARR = y_PARR[sel_PARR];d_PARR = d_PARR[sel_PARR]

                # Here I proceed only if the date is inside the Date_Num_limit fixed
                if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

                    # Here I extract the PARR along the isopycnal
                    sel_layer_PARR = (y_PARR >= reference_isopycnal_down) & (y_PARR < reference_isopycnal_up)
                    if np.sum(sel_layer_PARR) > 0:  # If sel_layer_PARR has some True values, then I take as the PARR of this isopycnal the mean of the PARR values in correspondence with these layers
                        PARR_tmp = np.mean(z_PARR[sel_layer_PARR])
                        depth_PARR_tmp2 = np.mean(d_PARR[sel_layer_PARR])
                        dens_PARR_tmp2 = np.mean(y_PARR[sel_layer_PARR])
                        PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                        depth_PARR_tmp = np.append(depth_PARR_tmp, depth_PARR_tmp2)
                        dens_PARR_tmp = np.append(dens_PARR_tmp, dens_PARR_tmp2)
                    else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the PARR at that depth
                        depth_PARR_tmp2 = np.array([])
                        dens_PARR_tmp2 = np.array([])
                        for iy in range(0, y_PARR.size - 1):
                            if y_PARR[iy] <= reference_isopycnal < y_PARR[iy + 1]:
                                dist = (reference_isopycnal - y_PARR[iy]) / (y_PARR[iy + 1] - y_PARR[iy])
                                PARR_tmp = z_PARR[iy] + (z_PARR[iy + 1] - z_PARR[iy]) * dist
                                d_tmp = d_PARR[iy] + (d_PARR[iy + 1] - d_PARR[iy]) * dist
                                PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                                depth_PARR_tmp2 = np.append(depth_PARR_tmp2, d_tmp)
                                dens_PARR_tmp2 = np.append(dens_PARR_tmp2, reference_isopycnal)
                        if depth_PARR_tmp2.size > 0:
                            depth_PARR_tmp = np.append(depth_PARR_tmp, np.mean(depth_PARR_tmp2))
                        if dens_PARR_tmp2.size > 0:
                            dens_PARR_tmp = np.append(dens_PARR_tmp, np.mean(dens_PARR_tmp2))

            PARR_list[i_isop] = np.mean(PARR_isopycnal)
            PARR_std_list[i_isop] = np.std(PARR_isopycnal)
            list_depth_PARR[i_isop] = np.mean(depth_PARR_tmp)
            list_dens_PARR[i_isop] = np.mean(dens_PARR_tmp)

        # I convert the PARR and the oxygen respiration rates (in micromolO2/kg/d) to the total amount of carbon consumption
        # between depth0 and depthf, and between day0 and dayf (in mgC/m2)
        # *Oxy2C -> to micromolC/kg/d
        # *mol2gC -> to microgC/kg/d
        # /1000 -> to mgC/kg/d
        # *density -> to mgC/m3/d
        # *layer_thickness*ndays -> to mgC/m2

        POC_resp_mgC_m2_vs_depth_list[:,iRespi] = PARR_list.copy()*list_dens_PARR*Oxy2C*mol2gC*layer_thickness*ndays/1000
        POC_resp_mgC_m2_std_vs_depth_list[:,iRespi] = PARR_std_list.copy()*list_dens_PARR*Oxy2C*mol2gC*layer_thickness*ndays/1000

    ########################################################################################################################
    # Here I calculate the carbon budget for depth0—depthf layer
    ########################################################################################################################
    Date_Num_Flux = x_filtered
    depth_POC_resp = list_depth_PARR
    depth_02_resp = depth_isopycnal

    ############### I calculate the integrated POC (MiP+MaP+bbp), between depth0 and depthf, for day0 and dayf. I transform it to mgC/m2

    # I extract the index of Integrated_POC_mgC_m3 which correspond to day0 (and dayf)
    tmp = list_dates_Integrated_POC - day0_float
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = list_dates_Integrated_POC - dayf_float
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]

    Integrated_POC_day0_mgC_m2 = Integrated_POC_mgC_m3[idx0] * layer_thickness
    Integrated_POC_dayf_mgC_m2 = Integrated_POC_mgC_m3[idxf] * layer_thickness
    Integrated_POC_day0_mgC_m2_std = Integrated_POC_mgC_m3_std[idx0] * layer_thickness
    Integrated_POC_dayf_mgC_m2_std = Integrated_POC_mgC_m3_std[idxf] * layer_thickness
    Delta_Integrated_POC = Integrated_POC_dayf_mgC_m2 - Integrated_POC_day0_mgC_m2
    Delta_Integrated_POC_std = np.sqrt( Integrated_POC_dayf_mgC_m2_std**2 + Integrated_POC_day0_mgC_m2_std**2 )

    Integrated_POC_extended_day0_mgC_m2 = Integrated_POC_extended_mgC_m3[idx0] * layer_thickness
    Integrated_POC_extended_dayf_mgC_m2 = Integrated_POC_extended_mgC_m3[idxf] * layer_thickness
    Integrated_POC_extended_day0_mgC_m2_std = Integrated_POC_extended_mgC_m3_std[idx0] * layer_thickness
    Integrated_POC_extended_dayf_mgC_m2_std = Integrated_POC_extended_mgC_m3_std[idxf] * layer_thickness
    Delta_Integrated_POC_extended = Integrated_POC_extended_dayf_mgC_m2 - Integrated_POC_extended_day0_mgC_m2
    Delta_Integrated_POC_extended_std = np.sqrt( Integrated_POC_extended_dayf_mgC_m2_std**2 + Integrated_POC_extended_day0_mgC_m2_std**2 )

    ############### I calculate the amount of POC entering from depht0 and exiting from dayf between day0 and dayf (in mgC/m2)

    # I extract the index of Flux_depth0/Flux_depthf which correspond to day0 (and dayf)
    tmp = Date_Num_Flux - day0_float
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = Date_Num_Flux - dayf_float
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]

    tmp = Flux_depth0[idx0:idxf]
    Flux_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_depthf[idx0:idxf]
    Flux_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_eta_b_depth0[idx0:idxf]
    Flux_eta_b_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_eta_b_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_eta_b_depthf[idx0:idxf]
    Flux_eta_b_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_eta_b_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_extended_depth0[idx0:idxf]
    Flux_extended_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_extended_depthf[idx0:idxf]
    Flux_extended_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_extended_eta_b_depth0[idx0:idxf]
    Flux_extended_eta_b_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_eta_b_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_extended_eta_b_depthf[idx0:idxf]
    Flux_extended_eta_b_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_eta_b_depthf_mgC_m2_std = np.std(tmp) * ndays

    Delta_flux = Flux_depth0_mgC_m2 - Flux_depthf_mgC_m2
    Delta_flux_eta_b = Flux_eta_b_depth0_mgC_m2 - Flux_eta_b_depthf_mgC_m2
    Delta_flux_extended = Flux_extended_depth0_mgC_m2 - Flux_extended_depthf_mgC_m2
    Delta_flux_extended_eta_b = Flux_extended_eta_b_depth0_mgC_m2 - Flux_extended_eta_b_depthf_mgC_m2

    Delta_flux_std = np.sqrt( Flux_depth0_mgC_m2_std**2 + Flux_depthf_mgC_m2_std**2 )
    Delta_flux_eta_b_std = np.sqrt( Flux_eta_b_depth0_mgC_m2_std**2 + Flux_eta_b_depthf_mgC_m2_std**2 )
    Delta_flux_extended_std = np.sqrt( Flux_extended_depth0_mgC_m2_std**2 + Flux_extended_depthf_mgC_m2_std**2 )
    Delta_flux_extended_eta_b_std = np.sqrt( Flux_extended_eta_b_depth0_mgC_m2_std**2 + Flux_extended_eta_b_depthf_mgC_m2_std**2 )

    Theoretical_Budget = Delta_flux - Delta_Integrated_POC
    Theoretical_Budget_eta_b = Delta_flux_eta_b - Delta_Integrated_POC
    Theoretical_Budget_extended = Delta_flux_extended - Delta_Integrated_POC_extended
    Theoretical_Budget_extended_eta_b = Delta_flux_extended_eta_b - Delta_Integrated_POC_extended

    Theoretical_Budget_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_eta_b_std = np.sqrt( Delta_flux_eta_b_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_extended_std**2 )
    Theoretical_Budget_extended_eta_b_std = np.sqrt( Delta_flux_extended_eta_b_std**2 + Delta_Integrated_POC_extended_std**2 )

    ############### I calculate the PARR and oxygen consumption between depth0 and depthf (in mgC/m2)

    # I extract the index of POC_resp_mgC_m2_vs_depth which correspond to depth0 (and depthf)
    tmp = depth_POC_resp - depth0
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = depth_POC_resp - depthf
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    POC_resp_mgC_m2_list=np.mean(POC_resp_mgC_m2_vs_depth_list[idx0:idxf,:],0)
    POC_resp_mgC_m2_std_list=np.mean(POC_resp_mgC_m2_std_vs_depth_list[idx0:idxf,:],0)

    # I extract the index of O2_resp_mgC_m2_vs_depth which correspond to depth0 (and depthf)
    tmp = depth_02_resp - depth0
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = depth_02_resp - depthf
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    O2_resp_mgC_m2=np.mean(O2_resp_mgC_m2_vs_depth[idx0:idxf])*-1
    O2_resp_mgC_m2_ci=np.mean(O2_resp_mgC_m2_vs_depth_ci[idx0:idxf,:],0)*-1
    O2_depth=np.mean(depth_02_resp[idx0:idxf])

    ############### I return the data
    return Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_eta_b,Theoretical_Budget_eta_b_std,\
           Theoretical_Budget_extended,Theoretical_Budget_extended_std,\
           Theoretical_Budget_extended_eta_b,Theoretical_Budget_extended_eta_b_std,\
           POC_resp_mgC_m2_list,POC_resp_mgC_m2_std_list,O2_resp_mgC_m2,O2_resp_mgC_m2_ci,O2_depth,list_Respi_types,n_profiles, \
           Delta_flux_eta_b, Delta_Integrated_POC,Delta_flux_eta_b_std, Delta_Integrated_POC_std

#######################################################################
# Parameters for the carbon budget calculation
#######################################################################
day0=datetime.datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime.datetime(2021,7,30)        # starting date for the carbon budget calculation
ndays=(dayf-day0).days          # number of days
dens00=1026.3                   # starting isopycnal
layer_thickness=0.1             # thickness of the layer considered (in kg/m3)
delta_dens=0.05                 # every time I do a loop, how much I do increase depth0
densff=1027.5                   # final isopycnal investigated

dens0_list=np.r_[dens00:densff-layer_thickness+0.1:delta_dens]

#######################################################################
# I loop on the different depths
#######################################################################
Theoretical_Budget_list_w1 = np.array([])
Theoretical_Budget_eta_b_list_w1 = np.array([])
Theoretical_Budget_extended_list_w1 = np.array([])
Theoretical_Budget_extended_eta_b_list_w1 = np.array([])
Theoretical_Budget_std_list_w1 = np.array([])
Theoretical_Budget_eta_b_std_list_w1 = np.array([])
Theoretical_Budget_extended_std_list_w1 = np.array([])
Theoretical_Budget_extended_eta_b_std_list_w1 = np.array([])
POC_resp_mgC_m2_list_w1 = np.array([])
POC_resp_mgC_m2_std_list_w1 = np.array([])
O2_resp_mgC_m2_list_w1 = np.array([])
O2_resp_mgC_m2_ci_list_w1 = np.array([])
O2_depth_list_w1 = np.array([])
dens0=dens0_list[0]
for dens0 in dens0_list:
    densf = dens0 + layer_thickness
    (Theoretical_Budget, Theoretical_Budget_std, Theoretical_Budget_eta_b, Theoretical_Budget_eta_b_std,
     Theoretical_Budget_extended, Theoretical_Budget_extended_std, Theoretical_Budget_extended_eta_b,
     Theoretical_Budget_extended_eta_b_std, POC_resp_mgC_m2, POC_resp_mgC_m2_std, O2_resp_mgC_m2, O2_resp_mgC_m2_ci,
     O2_depth, RespirationTypes, n_profiles, Delta_flux_eta_b, Delta_Integrated_POC, Delta_flux_eta_b_std,
     Delta_Integrated_POC_std) = carbon_budget_calculation(dens0, densf, day0, dayf)

    Theoretical_Budget_list_w1=np.append(Theoretical_Budget_list_w1,Theoretical_Budget)
    Theoretical_Budget_eta_b_list_w1=np.append(Theoretical_Budget_eta_b_list_w1,Theoretical_Budget_eta_b)
    Theoretical_Budget_extended_list_w1=np.append(Theoretical_Budget_extended_list_w1,Theoretical_Budget_extended)
    Theoretical_Budget_extended_eta_b_list_w1=np.append(Theoretical_Budget_extended_eta_b_list_w1,Theoretical_Budget_extended_eta_b)
    Theoretical_Budget_std_list_w1=np.append(Theoretical_Budget_std_list_w1,Theoretical_Budget_std)
    Theoretical_Budget_eta_b_std_list_w1=np.append(Theoretical_Budget_eta_b_std_list_w1,Theoretical_Budget_eta_b_std)
    Theoretical_Budget_extended_std_list_w1=np.append(Theoretical_Budget_extended_std_list_w1,Theoretical_Budget_extended_std)
    Theoretical_Budget_extended_eta_b_std_list_w1=np.append(Theoretical_Budget_extended_eta_b_std_list_w1,Theoretical_Budget_extended_eta_b_std)
    POC_resp_mgC_m2_list_w1=np.append(POC_resp_mgC_m2_list_w1,POC_resp_mgC_m2,axis=0)
    POC_resp_mgC_m2_std_list_w1=np.append(POC_resp_mgC_m2_std_list_w1,POC_resp_mgC_m2_std,axis=0)
    O2_resp_mgC_m2_list_w1=np.append(O2_resp_mgC_m2_list_w1,O2_resp_mgC_m2)
    O2_resp_mgC_m2_ci_list_w1=np.append(O2_resp_mgC_m2_ci_list_w1,O2_resp_mgC_m2_ci,axis=0)
    O2_depth_list_w1=np.append(O2_depth_list_w1,O2_depth)

O2_resp_mgC_m2_ci_list_w1=O2_resp_mgC_m2_ci_list_w1.reshape(dens0_list.size,2)
POC_resp_mgC_m2_list_w1=POC_resp_mgC_m2_list_w1.reshape(dens0_list.size,len(RespirationTypes))
POC_resp_mgC_m2_std_list_w1=POC_resp_mgC_m2_std_list_w1.reshape(dens0_list.size,len(RespirationTypes))

########################################################################################################################
######### Fig. 04a
########################################################################################################################
fs=10
width, height = 0.78, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.15, width, height])
plt.plot(O2_resp_mgC_m2_list_w1,O2_depth_list_w1, 'k')
plt.scatter(O2_resp_mgC_m2_list_w1,O2_depth_list_w1, c='black',s=5)
plt.fill_betweenx(O2_depth_list_w1, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
for iResp in range(2,3):
    plt.plot(POC_resp_mgC_m2_list_w1[:,iResp], dens0_list + layer_thickness / 2, c='b')

plt.fill_betweenx(dens0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                  POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m2_list_w1[:, 0], dens0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m2_list_w1[:, 3], dens0_list + layer_thickness / 2, c='b',linestyle='dashed')
plt.plot(POC_resp_mgC_m2_list_w1[:, 4], dens0_list + layer_thickness / 2, c='b',linestyle='dashed')
plt.plot(POC_resp_mgC_m2_list_w1[:, 5], dens0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
plt.plot(POC_resp_mgC_m2_list_w1[:, 6], dens0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
plt.plot(Theoretical_Budget_list_w1, dens0_list + layer_thickness / 2, c='red')
plt.scatter(Theoretical_Budget_list_w1, dens0_list + layer_thickness / 2, c='red', s=5)
plt.fill_betweenx(dens0_list + layer_thickness / 2, Theoretical_Budget_list_w1 - Theoretical_Budget_std_list_w1*0.5, Theoretical_Budget_list_w1 + Theoretical_Budget_std_list_w1*0.5,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nresp. rate')

plt.xlim(-570,7500)
plt.ylabel('Depth (m)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
ax.text(-0.05, 1.045, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v02/Fig04a_v02.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Fig. 04b
########################################################################################################################
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.15, width, height])
plt.plot(O2_resp_mgC_m2_list_w1,O2_depth_list_w1, 'k')
plt.scatter(O2_resp_mgC_m2_list_w1,O2_depth_list_w1, c='black',s=5)
plt.fill_betweenx(O2_depth_list_w1, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
for iResp in range(9,10):
    plt.plot(POC_resp_mgC_m2_list_w1[:,iResp], dens0_list + layer_thickness / 2, c='b')

plt.fill_betweenx(dens0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                  POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m2_list_w1[:, 7], dens0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m2_list_w1[:, 10], dens0_list + layer_thickness / 2, c='b',linestyle='dashed')
plt.plot(POC_resp_mgC_m2_list_w1[:, 11], dens0_list + layer_thickness / 2, c='b',linestyle='dashed')
plt.plot(POC_resp_mgC_m2_list_w1[:, 12], dens0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
plt.plot(POC_resp_mgC_m2_list_w1[:, 13], dens0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
plt.plot(Theoretical_Budget_extended_list_w1, dens0_list + layer_thickness / 2, c='red')
plt.scatter(Theoretical_Budget_extended_list_w1, dens0_list + layer_thickness / 2, c='red', s=5)
plt.fill_betweenx(dens0_list + layer_thickness / 2, Theoretical_Budget_extended_list_w1 - Theoretical_Budget_extended_std_list_w1*0.5, Theoretical_Budget_extended_list_w1 + Theoretical_Budget_extended_std_list_w1*0.5,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nresp. rate')

plt.xlim(-570,7500)
plt.ylabel('Depth (m)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
ax.text(-0.05, 1.045, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v02/Fig04b_v02.pdf' ,dpi=200)
plt.close()
# endregion










