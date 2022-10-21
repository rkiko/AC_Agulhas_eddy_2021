
########################################################################################################################
########################################################################################################################
########################################################################################################################
######### SUPPLEMENTARY FIG TEMPERATURE ANOMALY
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Supplementary Fig. Temperature anomaly
import scipy.io as sp
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
# I load the reference temperature
#######################################################################
mat_contents = sp.loadmat("%s/MatHydroProperties.mat" % (storedir))
# CT_Float=mat_contents["CT_Float"]
CT_Ref=mat_contents["CT_Ref"]

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

CT_Ref=CT_Ref.T[0:Date_Num.size,:]
# CT_Float=CT_Float.T[0:Date_Num.size,:]

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
cons_temp_anom = cons_temp - CT_Ref

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
cons_temp_anom=cons_temp_anom[sel_insideEddy,:]

#######################################################################
# I calculate the mixed layer depth
#######################################################################
from oceanpy import mixed_layer_depth
mld=np.array([])
i=0
for i in range(0,cons_temp.shape[0]):
    depth_tmp=depth[i,:]
    temp_tmp=cons_temp[i,:]
    # I exclude nan values
    sel_non_nan=(depth_tmp!=99999)&(temp_tmp!=99999)
    temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
    mld_tmp,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
    mld=np.append(mld,mld_tmp)

#######################################################################
# I plot
#######################################################################
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)

parameter_ylabel_list=['Temperature anom. ($^{\circ}$C)']
parameter_panellabel_list=[' ']
parameter_shortname_list=['cons_temp_anom']
ipar=0
if ipar==0:   parameter=cons_temp_anom.copy()

#I filter the profiles
parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([])
i=0
for i in range(0,parameter.shape[0]):
    z = parameter[i, :]
    sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999) & (~np.isnan(z))
    if sum(sel) > 0:
        z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
        z = savgol_filter(z, 5, 1)
        parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
        depth_parameter = np.concatenate((depth_parameter, y1))
        dens_parameter = np.concatenate((dens_parameter, y2))

# parameter_filtered[parameter_filtered<0]=0
dens_parameter[dens_parameter<1026]=1026
dens_parameter[dens_parameter>1027.5]=1027.5
# I define the x and y arrays for the contourf plot
x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),200)
y2_parameter = np.linspace(dens_parameter.min(), dens_parameter.max(), 200)
# I interpolate
x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
dens_interp_depth = griddata((Date_Num_parameter,depth_parameter), dens_parameter, (x_parameter_g, y_parameter_g), method="nearest")
x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y2_parameter)
parameter_interp_dens = griddata((Date_Num_parameter,dens_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")

#cons temp anom in the ML
mld_int = np.interp(x_parameter, Date_Num, mld)
parameter_mld = np.zeros((mld_int.size,))
i = 0
for i in range(0, mld_int.size):
    tmp = parameter_interp_depth[:, i]
    sel_mld = y1_parameter <= (mld_int[i] - 20)
    parameter_mld[i] = np.mean(tmp[sel_mld])

cons_temp_anom_ML = parameter_mld.copy()

# Parameter in the eddy core
dens0 = 1026.82
dens1 = 1027.2397618090454  # calculated at step 4 of Fig. 3a
parameter_eddy_core = np.zeros((mld_int.size,))
i = 0
for i in range(0, mld_int.size):
    tmp = parameter_interp_dens[:, i]
    sel_tmp = (y2_parameter >= dens0) & (y2_parameter < dens1)
    parameter_eddy_core[i] = np.mean(tmp[sel_tmp])

parameter_eddy_core = savgol_filter(parameter_eddy_core, 5, 1)
cons_temp_anom_eddy_core = parameter_eddy_core.copy()

# Parameter between ML and lower boundary of the eddy core
parameter_ML_eddy_core_down = np.zeros((mld_int.size,))
i = 0
for i in range(0, mld_int.size):
    tmp = parameter_interp_dens[:, i]
    sel_tmp = (y1_parameter > mld_int[i] ) & (y2_parameter < dens1)
    parameter_ML_eddy_core_down[i] = np.mean(tmp[sel_tmp])

parameter_ML_eddy_core_down = savgol_filter(parameter_ML_eddy_core_down, 5, 1)
cons_temp_anom_ML_eddy_core_down = parameter_ML_eddy_core_down.copy()

from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
argument = 'temp_anom_ML_0413_0923'
arg_value=np.mean(cons_temp_anom_ML)
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'temp_anom_eddy_core_0413_0923'
arg_value=np.mean(cons_temp_anom_eddy_core)
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'temp_anom_ML_eddy_core_down_0413_0923'
arg_value=np.mean(cons_temp_anom_ML_eddy_core_down)
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'temp_anom_minimum_0413_0923'
arg_value=np.mean(np.min(parameter_interp_depth,axis=0))
write_latex_data(filename,argument,'%0.2f' % arg_value)

########################################################
####### (Useless part for the moment) I plot: versus depth
########################################################

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = y1_parameter.min(),800
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, parameter_interp_depth,cmap='RdBu_r')
plt.plot(Date_Num,mld,'w')
# plot3 = ax.contour(x_parameter, y1_parameter, dens_interp_depth, levels=[1026.35],colors='white', linestyles='dashed', linewidths=1, zorder=10 )#,cmap='RdBu')
# fmt = {}
# strs = ['1026.35 kg/m$^3$']
# for l, s in zip(plot3.levels, strs):
#     fmt[l] = s
# ax.clabel(plot3, plot3.levels[::], inline=True, fmt=fmt, fontsize=10)

plt.gca().invert_yaxis()
plt.vlines(day_start_eddy_merging,ymin=0,ymax=700,color='w',linestyles='dashed')
plt.vlines(day_end_eddy_merging,ymin=0,ymax=700,color='w',linestyles='dashed')
# draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
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
# I don't save cos I prefer the plot by Artemis
# plt.savefig('../Plots/Supplementary/Anom_ConsTemp_v05.pdf',dpi=200)
# plt.close()

# endregion


########################################################################################################################
########################################################################################################################
########################################################################################################################
######### SUPPLEMENTARY FIG SATELLITE CHL TIME SERIES
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Supplementary Fig. Satellite Chl Time Series
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
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable


a_file = open("%s/an69/data_an69.pkl" % storedir, "rb")
data_an69 = pickle.load(a_file)
chl_inside_mean1=data_an69['chl_inside_mean1']
chl_outside_mean1=data_an69['chl_outside_mean1']
Date_Num_Eddy1=data_an69['Date_Num_Eddy1']
chl_inside_mean2=data_an69['chl_inside_mean2']
chl_outside_mean2=data_an69['chl_outside_mean2']
Date_Num_Eddy2=data_an69['Date_Num_Eddy2']
a_file.close()
# I define the end of the time series
day0_insideEddy=np.array([2021,4,13])
dayf_insideEddy=np.array([2021,9,23])
day0_insideEddy=matlab_datenum(day0_insideEddy)
dayf_insideEddy=matlab_datenum(dayf_insideEddy)

width, height = 0.8, 0.5
set_xlim_lower, set_xlim_upper = min(Date_Num_Eddy1.min(),Date_Num_Eddy2.min()), Date_Num_Eddy1.max()+10
set_ylim_lower, set_ylim_upper = 0,1.2#anom_meanE_meanOut.min()*1.1,chl_inside_mean.max()*1.1
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(set_xlim_lower, set_xlim_upper))
plt.plot(Date_Num_Eddy1,chl_inside_mean1,'tomato',label='Chl inside eddy 1')
plt.plot(Date_Num_Eddy1,chl_outside_mean1,'seagreen',label='Chl outside eddy 1')
plt.plot(Date_Num_Eddy2,chl_inside_mean2,'lightsalmon',label='Chl inside eddy 2')
plt.plot(Date_Num_Eddy2,chl_outside_mean2,'lime',label='Chl outside eddy 2')
plt.vlines([day0_insideEddy,dayf_insideEddy],set_ylim_lower, set_ylim_upper,colors='black',label='Float inside eddy',linewidth=3,linestyles='dashed')
plt.hlines([set_ylim_lower+0.025,set_ylim_upper-0.025],day0_insideEddy, dayf_insideEddy,colors='black',linewidth=3,linestyles='dashed')
# I set xticks
nxticks = 10
xticks = np.linspace(set_xlim_lower, set_xlim_upper, nxticks)
xticklabels = []
for i in xticks:
    tmp=matlab_datevec(i).astype(int)
    xticklabels.append(datetime.datetime(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=12)
plt.ylabel('Chlorophyll (mg/m$^3$)', fontsize=15)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/Fig_Main_v05/Supplementary/Supplementary_v05.pdf' ,dpi=200)
plt.close()

#######################################################################
# I save some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=0;print(matlab_datevec(Date_Num_Eddy1[i]).astype(int))
argument = 'chl_Eddy1_20201018'
arg_value=chl_inside_mean1[i]
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=73;print(matlab_datevec(Date_Num_Eddy1[i]).astype(int))
i=285;print(matlab_datevec(Date_Num_Eddy1[i]).astype(int))
argument = 'chl_Eddy1_0101to0108'
arg_value=np.mean(chl_inside_mean1[73:i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=175;print(matlab_datevec(Date_Num_Eddy1[i]).astype(int))
i=339;print(matlab_datevec(Date_Num_Eddy1[i]).astype(int))
argument = 'chl_percentage_difference_inside_outside'
arg_value = np.mean(chl_inside_mean1[175:i]-chl_outside_mean1[175:i])/np.mean(chl_outside_mean1[175:i])*100
write_latex_data(filename,argument,'%d' % np.round(arg_value))
i=3;print(matlab_datevec(Date_Num_Eddy2[i]).astype(int))
argument = 'chl_Eddy2_20200921'
arg_value=chl_inside_mean2[i]
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=105;print(matlab_datevec(Date_Num_Eddy2[i]).astype(int))
i=317;print(matlab_datevec(Date_Num_Eddy2[i]).astype(int))
argument = 'chl_Eddy2_0101to0108'
arg_value=np.mean(chl_inside_mean2[105:i])
write_latex_data(filename,argument,'%0.2f' % arg_value)

# endregion


########################################################################################################################
########################################################################################################################
########################################################################################################################
######### SUPPLEMENTARY FIG TIME SERIES OF TEMP, CHL, MiP etc. in ML and EDDY CORE
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Supplementary Fig. Part 1: Time Series of temp, chl, doxy, and bbp. in ML and eddy core
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
from oceanpy import bbp700toPOC_Koestner
sPOC=bbp700.copy()*0+99999
sPOC_Koestner=bbp700.copy()*0+99999
i=0
for i in range(0,bbp700.shape[0]):
    bbp700tmp=bbp700[i,:]
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
    chl_tmp=chla[i,:]
    # I exclude nan values
    sel=(bbp700tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)&(chl_tmp!=99999)
    bbp700tmp=bbp700tmp[sel]
    depth_tmp=depth_tmp[sel]
    temp_tmp=temp_tmp[sel]
    chl_tmp=chl_tmp[sel]
    # I convert to small POC (sPOC) and I set to 0 values <0
    sPOC_tmp = bbp700toPOC(bbp700tmp, depth_tmp, temp_tmp)
    sPOC_tmp[sPOC_tmp<0]=0
    sPOC[i,sel]=sPOC_tmp
    sPOC_tmp = bbp700toPOC_Koestner(bbp700tmp, chl_tmp)
    sPOC_tmp[np.isnan(sPOC_tmp)]=0
    sPOC_Koestner[i,sel]=sPOC_tmp

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
sPOC_Koestner=sPOC_Koestner[sel_insideEddy,:]

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
# I plot
#######################################################################
day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)

parameter_ylabel_list=['Temperature ($^{\circ}$C)','Chlorophyll $a$ (mg/m$^3$)','Dissolved oxygen ($\mu$mol/kg)','$b_{bp}$POC (mgC $m^{-3}$)','$b_{bp}$POC (mgC $m^{-3}$)']
parameter_panellabel_list=['a','c','b','e','e']
parameter_panellabel_list_EC=['a','c','b','c','c']
parameter_shortname_list=['cons_temp','chla','doxy','bbpPOC_Cetinic','bbpPOC_Koestner']
ipar=0
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0:   parameter=cons_temp.copy()
    elif ipar == 1: parameter=chla.copy()
    elif ipar == 2: parameter=doxy.copy()
    elif ipar == 3: parameter=sPOC.copy()
    elif ipar == 4: parameter=sPOC_Koestner.copy()

    #I filter the profiles
    parameter_filtered=np.array([]);Date_Num_parameter=np.array([]);depth_parameter=np.array([]);dens_parameter=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        z = parameter[i, :]
        sel=(z!=99999) & (depth[i,:]!=99999) & (dens[i,:]!=99999)
        if ipar == 3: sel = (sel) & (z <= 100)
        if sum(sel) > 0:
            z=z[sel];x=np.ones(z.shape)*Date_Num[i];y1=depth[i,sel];y2=dens[i,sel];y3=pres[i,sel]
            z = savgol_filter(z, 5, 1)
            parameter_filtered = np.concatenate((parameter_filtered, z));Date_Num_parameter = np.concatenate((Date_Num_parameter, x))
            depth_parameter = np.concatenate((depth_parameter, y1));dens_parameter = np.concatenate((dens_parameter, y2))

    parameter_filtered[parameter_filtered<0]=0
    # I define the x and y arrays for the contourf plot
    x_parameter = np.linspace(Date_Num_parameter.min(),Date_Num_parameter.max(),100)
    y1_parameter = np.linspace(depth_parameter.min(),depth_parameter.max(),200)
    y2_parameter = np.linspace(dens_parameter.min(), dens_parameter.max(), 200)
    # I interpolate
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y1_parameter)
    parameter_interp_depth = griddata((Date_Num_parameter,depth_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")
    x_parameter_g,y_parameter_g=np.meshgrid(x_parameter,y2_parameter)
    parameter_interp_dens = griddata((Date_Num_parameter,dens_parameter), parameter_filtered, (x_parameter_g, y_parameter_g), method="nearest")

    # Parameter in the mixed layer
    mld_int = np.interp(x_parameter, Date_Num, mld)
    parameter_mld = np.zeros((mld_int.size,))
    i = 0
    for i in range(0, mld_int.size):
        tmp = parameter_interp_depth[:, i]
        sel_mld = y1_parameter <= (mld_int[i] - 20)
        parameter_mld[i] = np.mean(tmp[sel_mld])

    if ipar == 0:
        temp_mld = parameter_mld.copy()
    elif ipar == 1:
        chla_mld = parameter_mld.copy()
    elif ipar == 2:
        doxy_mld = parameter_mld.copy()
    elif ipar == 3:
        bbpPOC_mld = parameter_mld.copy()
    elif ipar == 4:
        bbpPOC_Koestner_mld = parameter_mld.copy()

    width, height = 0.8, 0.7
    fig = plt.figure(1, figsize=(12, 4))
    ax = fig.add_axes([0.12, 0.4, width,height - 0.15])  # ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
    plt.plot(x_parameter, parameter_mld)
    plt.ylabel(parameter_ylabel_list[ipar])
    plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
    plt.xlim(x_parameter[0], x_parameter[-1])
    plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
    plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
    # I set xticks
    nxticks = 10
    xticks = np.linspace(Date_Num.min(), Date_Num.max(), nxticks)
    xticklabels = []
    for i in xticks:
        date_time_obj = date_reference + datetime.timedelta(days=i)
        xticklabels.append(date_time_obj.strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90, fontsize=14)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    # I save
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/TimeSeries_ML/TimeSeries_ML_%s_v05.pdf' % parameter_shortname_list[ipar], dpi=200)
    plt.close()

    if ipar == 1:
        chla_mld_integrated = parameter_mld.copy()*mld_int
        width, height = 0.8, 0.7
        fig = plt.figure(1, figsize=(12, 4))
        ax = fig.add_axes([0.12, 0.4, width,
                           height - 0.15])  # ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
        plt.plot(x_parameter, chla_mld_integrated)
        plt.ylabel('Integrated chlorophyll $a$ (mg/m$^2$)')
        plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
        plt.xlim(x_parameter[0], x_parameter[-1])
        plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
        plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
        # I set xticks
        nxticks = 10
        xticks = np.linspace(Date_Num.min(), Date_Num.max(), nxticks)
        xticklabels = []
        for i in xticks:
            date_time_obj = date_reference + datetime.timedelta(days=i)
            xticklabels.append(date_time_obj.strftime('%d %B'))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        plt.xticks(rotation=90, fontsize=14)
        # I add the panel label
        ax.text(-0.05, 1.05, 'd', transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        # I save
        plt.savefig('../Plots/Fig_Main_v05/Supplementary/TimeSeries_ML/TimeSeries_ML_chla_integrated_v05.pdf', dpi=200)
        plt.close()

        continue # I do not plot chl in the eddy core

    # Parameter in the eddy core
    dens0 = 1026.82
    dens1 = 1027.2397618090454  # calculated at step 4 of Fig. 3a
    parameter_eddy_core = np.zeros((mld_int.size,))
    i = 0
    for i in range(0, mld_int.size):
        tmp = parameter_interp_dens[:, i]
        sel_tmp = (y2_parameter >= dens0) & (y2_parameter < dens1)
        parameter_eddy_core[i] = np.mean(tmp[sel_tmp])

    parameter_eddy_core = savgol_filter(parameter_eddy_core, 5, 1)

    if ipar == 0:
        temp_eddy_core = parameter_eddy_core.copy()
    elif ipar == 2:
        doxy_eddy_core = parameter_eddy_core.copy()
    elif ipar == 3:
        bbpPOC_eddy_core = parameter_eddy_core.copy()
    elif ipar == 4:
        bbpPOC_Koestner_eddy_core = parameter_eddy_core.copy()

    fig = plt.figure(1, figsize=(12, 4))
    ax = fig.add_axes([0.12, 0.4, width, height - 0.15])
    plt.plot(x_parameter, parameter_eddy_core)
    plt.ylabel(parameter_ylabel_list[ipar])
    plt.xlim(x_parameter[0], x_parameter[-1])
    plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
    plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
    plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
    # I set xticks
    nxticks = 10
    xticks = np.linspace(Date_Num.min(), Date_Num.max(), nxticks)
    xticklabels = []
    for i in xticks:
        date_time_obj = date_reference + datetime.timedelta(days=i)
        xticklabels.append(date_time_obj.strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90, fontsize=14)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list_EC[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    # I save
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/TimeSeries_EddyCore/TimeSeries_EddyCore_%s_v05.pdf' % parameter_shortname_list[ipar], dpi=200)
    plt.close()

#######################################################################
# I save the some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=5;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'doxy_0421_eddy_core'
arg_value=np.mean(doxy_eddy_core[i])
write_latex_data(filename,argument,'%d' % arg_value)
i=74;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'doxy_0813_eddy_core'
arg_value=np.mean(doxy_eddy_core[i])
write_latex_data(filename,argument,'%d' % arg_value)
i=92;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'doxy_0912_eddy_core'
arg_value=np.mean(doxy_eddy_core[i])
write_latex_data(filename,argument,'%d' % arg_value)
i=10;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'bbp_April_eddy_core'
arg_value=np.mean(bbpPOC_eddy_core[0:i+1])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=48;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'bbp_0701to0923_eddy_core'
arg_value=np.mean(bbpPOC_eddy_core[i:])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=73;print(matlab_datevec(x_parameter[i]+matlab_datenum(1950,1,1)).astype(int))
argument = 'bbp_0812_eddy_core'
arg_value=np.mean(bbpPOC_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)



# endregion
# region Supplementary Fig. Part 2: Time Series of MiP MaP in ML and eddy core
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
dens=np.array(data_ecopart['Potential density [kg/m3]'][sel_filename])

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
parameter_shortname_list=['MiP_POC','MaP_POC']
parameter_panellabel_list=['f','g']
parameter_panellabel_list_EC=['d','e'] # Panel label for the eddy core
parameter_ylabel_list=['MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)']
max_parameter_list=np.array([2.15,0.30])
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar == 0: parameter=MiP_POC.copy()
    elif ipar == 1: parameter=MaP_POC.copy()

    parameter_filtered=np.array([]);Date_Num_filtered=np.array([]);depth_filtered=np.array([]);dens_filtered=np.array([])
    # I filter the prophiles
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel];k=dens[sel]
        z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2];k2=k[sel2]
        if sum(sel2)>0:
            z=savgol_filter(z,5,1)
            parameter_filtered = np.concatenate((parameter_filtered, z))
            Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
            depth_filtered = np.concatenate((depth_filtered, y2))
            dens_filtered = np.concatenate((dens_filtered, k2))

    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y1_parameter = np.linspace(depth_filtered.min(), depth_filtered.max(), 200)
    y2_parameter = np.linspace(dens_filtered.min(), dens_filtered.max(), 200)

    # I interpolate
    x_parameter_g, y_parameter_g = np.meshgrid(x_filtered, y1_parameter)
    parameter_interp_depth = griddata((Date_Num_filtered, depth_filtered), parameter_filtered,(x_parameter_g, y_parameter_g), method="nearest")
    x_parameter_g, y_parameter_g = np.meshgrid(x_filtered, y2_parameter)
    parameter_interp_dens = griddata((Date_Num_filtered, dens_filtered), parameter_filtered,(x_parameter_g, y_parameter_g), method="nearest")

    # I load the mixed layer
    a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an68/data_MLD_an68.pkl" % (home), "rb")
    data_an68 = pickle.load(a_file)
    mld = data_an68['mld']
    mld_date = data_an68['Date_Num']
    a_file.close()
    mld = mld[sel_insideEddy]
    mld_date = mld_date[sel_insideEddy]
    # I transform the x_filtered in matlab datenum format
    x_filtered_datenum = x_filtered.copy()
    i =0
    for i in range(0,x_filtered.size):
        d = datetime.datetime.utcfromtimestamp(x_filtered[i])
        x_filtered_datenum[i] = matlab_datenum(d.year,d.month,d.day,d.hour,d.minute,d.second)
        # datetime.utcfromtimestamp(Date_Num[i])

    # I interpolate the mixed layer
    mld_int = np.interp(x_filtered_datenum, mld_date, mld)

    # Parameter in the mixed layer
    parameter_mld = np.zeros((mld_int.size,))
    i = 0
    for i in range(0, mld_int.size):
        tmp = parameter_interp_depth[:, i]
        sel_mld = y1_parameter <= (mld_int[i] - 20)
        parameter_mld[i] = np.mean(tmp[sel_mld])

    if ipar == 0:
        MiP_mld = parameter_mld.copy()
    elif ipar == 1:
        MaP_mld = parameter_mld.copy()

    width, height = 0.8, 0.7
    fig = plt.figure(1, figsize=(12, 4))
    ax = fig.add_axes([0.12, 0.4, width,height - 0.15])  # ylim=(set_ylim_lower, set_ylim_upper),xlim=(Date_Num.min(), Date_Num.max()))
    plt.plot(x_filtered, parameter_mld)
    plt.ylabel(parameter_ylabel_list[ipar])
    plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
    plt.xlim(x_filtered[0], x_filtered[-1])
    plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
    plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[1], ymax=ax.get_ylim()[0], color='k')
    # I set xticks
    nxticks = 10
    xticks = np.linspace(x_filtered.min(), x_filtered.max(), nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90, fontsize=14)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    # I save
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/TimeSeries_ML/TimeSeries_ML_%s_v05.pdf' % parameter_shortname_list[ipar], dpi=200)
    plt.close()

    # Parameter in the eddy core
    dens0 = 1026.82
    dens1 = 1027.2397618090454  # calculated at step 4 of Fig. 3a
    parameter_eddy_core = np.zeros((mld_int.size,))
    i = 0
    for i in range(0, mld_int.size):
        tmp = parameter_interp_dens[:, i]
        sel_tmp = (y2_parameter >= dens0) & (y2_parameter < dens1)
        parameter_eddy_core[i] = np.mean(tmp[sel_tmp])

    parameter_eddy_core = savgol_filter(parameter_eddy_core, 5, 1)

    if ipar == 0:
        MiP_eddy_core = parameter_eddy_core.copy()
    elif ipar == 1:
        MaP_eddy_core = parameter_eddy_core.copy()

    fig = plt.figure(1, figsize=(12, 4))
    ax = fig.add_axes([0.12, 0.4, width, height - 0.15])
    plt.plot(x_filtered, parameter_eddy_core)
    plt.ylabel(parameter_ylabel_list[ipar])
    plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1])
    plt.xlim(x_filtered[0], x_filtered[-1])
    plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
    plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
    # I set xticks
    nxticks = 10
    xticks = np.linspace(x_filtered.min(), x_filtered.max(), nxticks)
    xticklabels = []
    for i in xticks:
        xticklabels.append(datetime.datetime.utcfromtimestamp(i).strftime('%d %B'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=90, fontsize=14)
    # I add the panel label
    ax.text(-0.05, 1.05, parameter_panellabel_list_EC[ipar], transform=ax.transAxes,fontsize=24, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
    # I add the grid
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    # I save
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/TimeSeries_EddyCore/TimeSeries_EddyCore_%s_v05.pdf' % parameter_shortname_list[ipar], dpi=200)
    plt.close()

#######################################################################
# I save the some values for the latex document
#######################################################################
from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
i=0;print(matlab_datevec(x_filtered_datenum[i]).astype(int))
argument = 'MiP_0413_eddy_core'
arg_value=np.mean(MiP_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=45;print(matlab_datevec(x_filtered_datenum[i]).astype(int))
argument = 'MiP_0626_eddy_core'
arg_value=np.mean(MiP_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
i=81;print(matlab_datevec(x_filtered_datenum[i]).astype(int))
argument = 'MiP_0825_eddy_core'
arg_value=np.mean(MiP_eddy_core[i])
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'MaP_0413to0923_eddy_core'
arg_value=np.mean(MaP_eddy_core)
write_latex_data(filename,argument,'%0.2f' % arg_value)
argument = 'MaP_0413to0923_eddy_core_std'
arg_value=np.std(MaP_eddy_core)
write_latex_data(filename,argument,'%0.2f' % arg_value)

# endregion


########################################################################################################################
########################################################################################################################
########################################################################################################################
######### SUPPLEMENTARY FIG PSD Extended + BBP POC
########################################################################################################################
########################################################################################################################
########################################################################################################################
#region Supplementary Fig. PSD Extended + BBP POC
import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from lin_fit import lin_fit

########################################################################################################################
#Time period that I plot
########################################################################################################################
day0=datetime(2021, 4,13)        # starting date
dayf=datetime(2021, 7,31)        # final date
day0_float = calendar.timegm(day0.timetuple())
dayf_float = calendar.timegm(dayf.timetuple())
ndays = (dayf - day0).days  # number of days

delta_bin=20 #thickness (in meters) of the bin I used to calculate the psd distribution
depth_f=600 #final depth (in meters) at which I calculate the psd distribution

list_depths_plot=np.r_[200:depth_f+delta_bin*0.5:delta_bin]

########################################################################################################################
#I process the data
########################################################################################################################
data = pd.read_csv("../Data/Ecopart_processed_data_356.tsv", sep='\t')

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
ic0 = data.columns.get_loc('Part conc frac [#/l] (ESD: 0.102-0.128 mm)')        # Starting index of the columns containing the size classes for which UVP works, and that I use for the fit
icf = data.columns.get_loc('Part conc frac [#/l] (ESD: 2.58-3.25 mm)')          # End index of the columns containing the size classes for which UVP works, and that I use for the fit
ice0 = data.columns.get_loc('Part conc frac [#/l] (ESD: 0.0254-0.032 mm)')      # Starting index of the columns containing the size classes for which UVP does not works, and for which I extrapolate the abundance
icef = ic0                                        # End index of the columns containing the size classes for which UVP does not works, and for which I extrapolate the abundance
PSD_columns=data.columns[ic0:icf]                 # Column labels containing the size classes for which UVP works, and that I use for the fit
PSD_columns_extrapolation=data.columns[ice0:icef] # Column labels containing the size classes for which UVP does not works, and for which I extrapolate the abundance

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:depth.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
filter_dates=(Date_Num>=day0_float)&(Date_Num<=dayf_float)

########################################################################################################################
#I calculate the size of the bins which I need to normalise the PSD distribution
########################################################################################################################

PSD_bin_start=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_end = np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_width=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_mean=np.squeeze(np.zeros((PSD_columns.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns.size):
    PSD_columns_name = PSD_columns[i_PSD]
    PSD_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_bin_width[i_PSD] = PSD_bin_end[i_PSD] - PSD_bin_start[i_PSD]
    PSD_bin_mean[i_PSD] = (PSD_bin_end[i_PSD] + PSD_bin_start[i_PSD])*0.5

########################################################################################################################
#I extract the width and mean values for the smaller size classes
########################################################################################################################

PSD_extrapolated_bin_start=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_end = np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_width=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_mean=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns_extrapolation.size):
    PSD_columns_name = PSD_columns_extrapolation[i_PSD]
    PSD_extrapolated_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_extrapolated_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_extrapolated_bin_width[i_PSD] = PSD_extrapolated_bin_end[i_PSD] - PSD_extrapolated_bin_start[i_PSD]
    PSD_extrapolated_bin_mean[i_PSD] = (PSD_extrapolated_bin_end[i_PSD] + PSD_extrapolated_bin_start[i_PSD])*0.5

PSD_extrapolated_bin_mean_log=np.log(PSD_extrapolated_bin_mean)

########################################################################################################################
#I do the loop for each profile, each depth, and each size class
########################################################################################################################
PSD_extrapolated=data.values[:,ice0:icef].copy()

ct_ok=0;ct_notOk=0
i_profile = 0
for i_profile in range (0,list_dates.size):
    sel_dates = Date_Num == list_dates[i_profile]
    list_depths = depth[sel_dates]
    i_depth = 1
    for i_depth in range(0, list_depths.size):
        PSD_depth = np.squeeze(np.zeros((PSD_columns.size,1)))
        i_PSD = 0
        for i_PSD in range(0, PSD_columns.size):
            PSD_columns_name = PSD_columns[i_PSD]
            PartConc = np.array(data[PSD_columns_name][sel_filename])  # I select a given size class data and I exclude the parkings
            PartConc = PartConc[sel_dates]                             # I select the data for the specific profile (i.e., date)
            PSD_depth[i_PSD] = PartConc[i_depth]                         # I select the data for the specific depth

        PSD_depth = PSD_depth/PSD_bin_width             # I divide by PSD_bin_width in order to normalise
        sel = (~np.isnan(PSD_depth)) & (PSD_depth != 0) # I exlcude nan and zero values
        PSD_depth=PSD_depth[sel]
        #I extrapolate the concentration only if the size distribution has at least 3 values
        if PSD_depth.size>2:
            PSD_bin_mean_tmp = PSD_bin_mean[sel]
            x = np.log(PSD_bin_mean_tmp)
            y = np.log(PSD_depth)
            sel=~np.isinf(y);x=x[sel];y=y[sel]              # I exclude inf values
            (interpol, slpe_ci, _, signif, signif_label) = lin_fit(x,y)
            #I extrapolate the concentration only if the fit is significant or the pvalue is less than 0.05
            if (signif>0)|(interpol.pvalue<=0.05):
                ct_ok = ct_ok +1
                PSD_extrapolated_depth_tmp = np.exp(PSD_extrapolated_bin_mean_log * interpol.slope + interpol.intercept)
                PSD_extrapolated[ np.squeeze(np.where(sel_filename))[sel_dates][i_depth] ,: ] = PSD_extrapolated_depth_tmp*PSD_extrapolated_bin_width
            else:
                ct_notOk = ct_notOk +1
        else:
            ct_notOk = ct_notOk + 1

PSD_extrapolated = PSD_extrapolated.astype('float32')

for i in range(0, PSD_columns_extrapolation.__len__()):
    data[PSD_columns_extrapolation[i]] = PSD_extrapolated[:, i]




########################################################################################################################
########################################################################################################################
########################################################################################################################
#I convert each size class concentration in POC and flux
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
depth=depth[filter_dates]
Date_Num=Date_Num[filter_dates]

########################################################################################################################
#I import the functions necessary
########################################################################################################################

from paruvpy import poc_flux_func_1_class
from paruvpy import poc_cont_func_1_class

########################################################################################################################
#I extract the width and mean values for the all the size classes
########################################################################################################################
PSD_columns_all=data.columns[ice0:icf]

PSD_all_bin_start=np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_end = np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_width=np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_mean=np.squeeze(np.zeros((PSD_columns_all.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns_all.size):
    PSD_columns_name = PSD_columns_all[i_PSD]
    PSD_all_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_all_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_all_bin_width[i_PSD] = PSD_all_bin_end[i_PSD] - PSD_all_bin_start[i_PSD]
    PSD_all_bin_mean[i_PSD] = (PSD_all_bin_end[i_PSD] + PSD_all_bin_start[i_PSD])*0.5

########################################################################################################################
# I calculate the POC and plot it. I add the bbp POC for comparison
########################################################################################################################
bbp=np.array(data['bbp POC [mgC/m3]'][sel_filename][filter_dates])
bbp_Koestner=np.array(data['bbp POC Koestner [mgC/m3]'][sel_filename][filter_dates])

i_depth=0
for i_depth in range(0,list_depths_plot.size-1):
    depth0_tmp = list_depths_plot[i_depth]
    depthf_tmp = list_depths_plot[i_depth + 1]
    depth_tmp = (depthf_tmp + depth0_tmp) * 0.5
    PartConc = data[PSD_columns_all][sel_filename][filter_dates]
    sel_depth = (depth >= depth0_tmp) & (depth < depthf_tmp)
    PartConc_depth = PartConc[sel_depth]
    PartConc_depth = np.mean(PartConc_depth,0)
    bbp_depth = np.mean(bbp [sel_depth])
    bbp_depth_std = np.std(bbp [sel_depth])
    bbp_Koestner_depth = np.mean(bbp_Koestner [sel_depth])
    bbp_Koestner_depth_std = np.std(bbp_Koestner [sel_depth])

    Flux=np.squeeze(np.zeros((1,PSD_columns_all.size)))
    POC=np.squeeze(np.zeros((1,PSD_columns_all.size)))
    i_PSD = 0
    for i_PSD in range(0,PSD_columns_all.size):
        PSD_columns_name = PSD_columns_all[i_PSD]
        Flux[i_PSD]= PartConc_depth[i_PSD] * poc_flux_func_1_class(PSD_columns_name)
        POC[i_PSD] = PartConc_depth[i_PSD] * poc_cont_func_1_class(PSD_columns_name)

    width, height = 0.75, 0.75
    set_ylim_lower, set_ylim_upper = 0.5*10**-3, 5*10 ** 3
    set_xlim_lower, set_xlim_upper = 0.0009, 2.7
    a=0
    b=5
    ##############
    #POC plot
    ##############
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.2, 0.15, width, height], xlim=(set_xlim_lower, set_xlim_upper), ylim=(set_ylim_lower, set_ylim_upper))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(PSD_all_bin_mean[b:-1], POC[b:-1], '--b.', label='UVP',linewidth=1)
    plt.plot(PSD_all_bin_mean[a:b+1], POC[a:b+1], c='b', label='UVP (extrap.)')
    # plt.scatter(0.013, bbp_depth, label='bbp Cetinic',c='red')
    # plt.errorbar(0.013, bbp_depth, xerr=0.012,yerr=bbp_depth_std, capsize=5,c='red')
    plt.scatter(0.013, bbp_Koestner_depth, label='b$_{bp}$POC',c='darkgoldenrod')
    plt.errorbar(0.013, bbp_Koestner_depth, xerr=0.012,yerr=bbp_depth_std, capsize=5,c='darkgoldenrod')
    plt.xlabel('Size (mm)', fontsize=10)
    plt.ylabel('POC (mgC/m$^3$)', fontsize=8)
    if i_depth==0:
        plt.legend(fontsize=10)
    plt.title('Depth: %d m' % depth_tmp, fontsize=10)#; R$^2$=%0.2f, Signif: %s' % (depth_tmp, interpol.rvalue ** 2, signif_label)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/ExtendedPSD/ExtendedPSD_%dm_v05.pdf' % (depth_tmp), dpi=200)
    plt.close()
    ##############
    #POC plot normalised
    ##############
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.2, 0.15, width, height], xlim=(set_xlim_lower, set_xlim_upper), ylim=(set_ylim_lower, set_ylim_upper))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(PSD_all_bin_mean[b:-1], POC[b:-1]/PSD_all_bin_mean[b:-1], '--b.', label='UVP',linewidth=1)
    plt.plot(PSD_all_bin_mean[a:b+1], POC[a:b+1]/PSD_all_bin_mean[a:b+1], c='b', label='UVP (extrap.)')
    # plt.scatter(0.013, bbp_depth/0.024, label='bbp Cetinic',c='red')
    # plt.errorbar(0.013, bbp_depth/0.024, xerr=0.012,yerr=bbp_depth_std, capsize=5,c='red')
    plt.scatter(0.013, bbp_Koestner_depth/0.024, label='b$_{bp}$POC',c='darkgoldenrod')
    plt.errorbar(0.013, bbp_Koestner_depth/0.024, xerr=0.012,yerr=bbp_depth_std, capsize=5,c='darkgoldenrod')
    plt.xlabel('Size (mm)', fontsize=10)
    plt.ylabel('Normalised POC (mgC/m$^3$/mm)', fontsize=8)
    if i_depth==0:
        plt.legend(fontsize=10)
    plt.title('Depth: %d m' % depth_tmp, fontsize=10)#; R$^2$=%0.2f, Signif: %s' % (depth_tmp, interpol.rvalue ** 2, signif_label)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/Fig_Main_v05/Supplementary/ExtendedPSD_Norm/ExtendedPSD_Norm_%dm_v05.pdf' % (depth_tmp), dpi=200)
    plt.close()




#endregion
