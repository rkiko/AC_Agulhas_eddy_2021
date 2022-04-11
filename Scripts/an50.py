import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from datetime import datetime
import calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename_mouw='%s/GIT/AC_Agulhas_eddy_2021/Data/an24/GO_flux.tab' % home

data=pd.read_csv(filename_mouw, sep='\t', header=87)
lon=data.Longitude
lat=data.Latitude
POC=data['POC flux [mg/m**2/day]']
Depth_trap=data['Depth water [m] (Sediment trap deployment depth)']
StartDate_trap=data['Date/Time (Deployed)']


########################################################################################################################
########################################################################################################################
# Flux from our BGC Argo float
########################################################################################################################
########################################################################################################################

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
day0=datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime(2021,9,24)       # final date for the carbon budget calculation
delta_depth=15                  # around of the depth which I consider when extracting the flux
day0_float=calendar.timegm(day0.timetuple())
dayf_float=calendar.timegm(dayf.timetuple())

########################################################################################################################
# I process the data
########################################################################################################################
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
Flux = np.array(data['Flux_mgC_m2'][sel_filename])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I select the data only in the prescribed period
list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]

##############################################
# I filter the Flux

ipar=0
for ipar in range(0,4):
    if ipar==0: parameter=Flux.copy()
    if ipar==1: parameter=Flux_eta_b.copy()
    if ipar==2: parameter=Flux_extended.copy()
    if ipar==3: parameter=Flux_extended_eta_b.copy()
    # I filter the flux prophiles
    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
        z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2)>0:
            z=savgol_filter(z,5,1)
            parameter_filtered = np.concatenate((parameter_filtered, z))
            Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
            depth_filtered = np.concatenate((depth_filtered, y2))

    # I define the x and y arrays for the interpolation
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    #I extract the flux
    x = x_filtered.copy()
    y = y_filtered.copy()
    z = parameter_interp.copy()

    depthf = 600  # final depth
    sel_depthf_600 = (np.abs(y) > depthf - delta_depth) & (np.abs(y) <= depthf + delta_depth)
    Flux_filtered_tmp = z[sel_depthf_600,:]
    if ipar==0: Flux_filtered_depthf_600 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==1: Flux_eta_b_filtered_depthf_600 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==2: Flux_extended_filtered_depthf_600 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==3: Flux_extended_eta_b_filtered_depthf_600 = np.mean(Flux_filtered_tmp,axis=0)

    depthf = 200  # final depth
    sel_depthf_200 = (np.abs(y) > depthf - delta_depth) & (np.abs(y) <= depthf + delta_depth)
    Flux_filtered_tmp = z[sel_depthf_200,:]
    if ipar==0: Flux_filtered_depthf_200 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==1: Flux_eta_b_filtered_depthf_200 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==2: Flux_extended_filtered_depthf_200 = np.mean(Flux_filtered_tmp,axis=0)
    if ipar==3: Flux_extended_eta_b_filtered_depthf_200 = np.mean(Flux_filtered_tmp,axis=0)


########################################################################################################################
########################################################################################################################
# I plot flux at 200 and 600 m
########################################################################################################################
########################################################################################################################
fs=9
width, height = 0.82, 0.8

# First plot: flux calculated without considering smallest size classes and with old eta and b values
fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x,Flux_filtered_depthf_200,'r',label='Flux at 200 m')
plt.plot(x,Flux_filtered_depthf_600,'b',label='Flux at 600 m')
plt.xlim(x.min(),x.max())
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
plt.savefig('../Plots/an50/01oldway_POCFlux_200_600m_%d%02d%02dto%d%02d%02d_an50.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Second plot: flux calculated without considering smallest size classes and with new eta and b values
fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x,Flux_eta_b_filtered_depthf_200,'r',label='Flux at 200 m')
plt.plot(x,Flux_eta_b_filtered_depthf_600,'b',label='Flux at 600 m')
plt.xlim(x.min(),x.max())
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
plt.savefig('../Plots/an50/02eta_b_POCFlux_200_600m_%d%02d%02dto%d%02d%02d_an50.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Third plot: flux calculated considering smallest size classes but with old eta and b values
fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x,Flux_extended_eta_b_filtered_depthf_200,'r',label='Flux at 200 m')
plt.plot(x,Flux_extended_eta_b_filtered_depthf_600,'b',label='Flux at 600 m')
plt.xlim(x.min(),x.max())
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
plt.savefig('../Plots/an50/03extended_POCFlux_200_600m_%d%02d%02dto%d%02d%02d_an50.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Fourth plot: flux calculated considering smallest size classes and with different eta and b values (based on the
# relationship between sinking speed and particle size)
fig = plt.figure(1, figsize=(5.5, 1.0))
ax = fig.add_axes([0.12, 0.1, width, height])
plt.plot(x,Flux_filtered_depthf_200,'r',label='Flux at 200 m')
plt.plot(x,Flux_filtered_depthf_600,'b',label='Flux at 600 m')
plt.xlim(x.min(),x.max())
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
plt.savefig('../Plots/an50/04extended_eta_b_200_600m_literature_%d%02d%02dto%d%02d%02d_an50.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()



