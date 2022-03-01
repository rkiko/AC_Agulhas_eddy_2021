import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import seawater as sw
import gsw
import pickle
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

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
PARR_nmol_l_h=np.array(data['Respi_nmolO2_l_h'][sel_filename])
temp=np.array(data['Temperature [degrees Celsius]'][sel_filename])
psal=np.array(data['Practical salinity [psu]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

#I tranform the pressure to depth
mask_depth=pressure!=99999 #I select only valid values
lat_tmp=lat[mask_depth]
pressure_tmp=pressure[mask_depth]
depth_tmp=sw.eos80.dpth(pressure_tmp, lat_tmp)
depth=np.ones(pressure.shape)*99999
depth[mask_depth]=depth_tmp

#I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
#salinity and temperature first
mask_dens=np.logical_and(pressure!=99999,temp!=99999,psal!=99999) # I exclude the points with value = 99999
lat_tmp=lat[mask_dens]
lon_tmp=lon[mask_dens]
pres_tmp=pressure[mask_dens]
psal_tmp=psal[mask_dens]
temp_tmp=temp[mask_dens]
abs_psal_tmp=gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp) # I compute absolute salinity
cons_tmp=gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)          # I compute conservative temperature
dens_tmp=gsw.density.sigma0(abs_psal_tmp, cons_tmp)
dens=np.ones(temp.shape)*99999
dens[mask_dens]=dens_tmp+1000

#I convert the PARR measured in micromol/kg/day
PARR_micromol_kg_day=-PARR_nmol_l_h/1000*24/(dens/1000)

list_dates=np.unique(Date_Num)
parameter=PARR_micromol_kg_day.copy()
ipar=0
parameter_ylabel_list=['PARR ($\mu$mol kg$^{-}$ $d^{-1}$)']
max_parameter_list=np.array([6,65,0.6,1.15,0.012])
#for ipar in range(0,parameter_ylabel_list.__len__()):

# I filter the profiles
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

# I define the x and y arrays for the contourf plot
x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
# I interpolate
parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

##########################################
# I import the oxygen respiration rates calculated in an10
##########################################
a_file = open("%s/an10/data_an10.pkl" % storedir, "rb")
data_an10 = pickle.load(a_file)
x_doxy_RR=data_an10['x_parameter']
y_doxy_RR=data_an10['y1_parameter']
doxy_RR=data_an10['doxy_RR_interp_depth']
##########################################
##########################################

################################################################################################################
####### I plot
################################################################################################################
# Parameters for the plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = 200,600
lim_min=-0.07
lim_max=0.0
nlevels=20

################# Plot part: PARR as contourfilled, doxy_RR as contour lines
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
parameter_plot=parameter_interp.copy()
parameter_plot[parameter_plot>lim_max]=lim_max
parameter_plot[parameter_plot<lim_min]=lim_min
doxy_RR_plot=doxy_RR.copy()
doxy_RR_plot[doxy_RR_plot>lim_max]=lim_max
doxy_RR_plot[doxy_RR_plot<lim_min]=lim_min
plot1 = plt.contourf(x_filtered, y_filtered, parameter_plot,levels=nlevels,cmap='Blues_r')
plot2 = ax.contour(x_filtered, y_doxy_RR, doxy_RR_plot,levels=5,colors='black',linestyles='solid',linewidths=1)#,cmap='Blues_r')
ax.clabel(plot2, inline=1, fontsize=10)
plt.ylim(set_ylim_lower, set_ylim_upper)
plt.gca().invert_yaxis()
# draw colorbar
cbar = plt.colorbar(plot1)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.title('PARR as contour filled, doxy RR as contour lines', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an11/PARR_vs_time_and_depth_an11.pdf',dpi=200)
plt.close()

################# Plot part: doxy_RR as contourfilled, PARR as contour lines
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
parameter_plot=parameter_interp.copy()
parameter_plot[parameter_plot>lim_max]=lim_max
parameter_plot[parameter_plot<lim_min]=lim_min
doxy_RR_plot=doxy_RR.copy()
doxy_RR_plot[doxy_RR_plot>lim_max]=lim_max
doxy_RR_plot[doxy_RR_plot<lim_min]=lim_min
plot1 = plt.contourf(x_filtered, y_doxy_RR, doxy_RR_plot,levels=nlevels,cmap='Blues_r')
plot2 = ax.contour(x_filtered, y_filtered, parameter_plot,levels=5,colors='black',linestyles='solid',linewidths=1)#,cmap='Blues_r')
ax.clabel(plot2, inline=1, fontsize=10)
plt.ylim(set_ylim_lower, set_ylim_upper)
plt.gca().invert_yaxis()
# draw colorbar
cbar = plt.colorbar(plot1)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.title('doxy RR as contour filled, PARR as contour lines', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an11/Doxy_RR_vs_time_and_depth_an11.pdf',dpi=200)
plt.close()



