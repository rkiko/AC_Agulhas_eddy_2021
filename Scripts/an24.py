import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import cartopy
import pickle
from datetime import datetime
import calendar
from scipy.signal import savgol_filter
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
#Parameters to identify the area for which we want to extract the POC flux data
########################################################################################################################
lon0,lon1 = -10, 20
lat0, lat1 = -40,-20
Depth_max_trap = 60000
########################################################################################################################
#I extract the data according to the domain selection
########################################################################################################################
sel_in_domain=(lon>=lon0)&(lon<=lon1)&(lat>=lat0)&(lat<=lat1)&(Depth_trap<=Depth_max_trap)
print('Found %d traps' % sum(sel_in_domain))

POC_in_domain=POC[sel_in_domain]
lon_in_domain=lon[sel_in_domain]
lat_in_domain=lat[sel_in_domain]
Depth_trap_in_domain=Depth_trap[sel_in_domain]
StartDate_trap_in_domain = StartDate_trap[sel_in_domain]

print('Mean Flux is %0.1f mgC/m2/d' % np.mean(POC_in_domain))
print('Mean sediment trap depth is %0.0f' % np.mean(Depth_trap_in_domain))

########################################################################################################################
#I plot the flux vs depth
########################################################################################################################
width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(12, 8))
ax = fig.add_axes([0.12, 0.2, width, height])
plot = plt.scatter(POC_in_domain, Depth_trap_in_domain, c='blue')
plt.xlabel('POC Flux (mgC/m$^2$/d)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.ylim(0,4000),plt.xlim(0,np.max(POC_in_domain)*1.1)
plt.gca().invert_yaxis()
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an24/POC_vs_Depth_an24.pdf', dpi=200)
plt.close()

########################################################################################################################
#I plot the flux trap positions
########################################################################################################################

delta_plot=1
set_xlim_lower, set_xlim_upper = lon0-delta_plot,lon1+delta_plot
set_ylim_lower, set_ylim_upper = lat0-delta_plot,lat1+delta_plot
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.plot(lon_in_domain, lat_in_domain, 'k.', markersize=12)
ax.plot([5,18,18,5,5],[-36,-36,-32,-32,-36],'k',label='BGC Argo region')
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
gl.xlabel_style = {'fontsize': 15}
gl.ylabel_style = {'fontsize': 15}
ax.set_xlabel('Longitude', fontsize=15)
ax.set_ylabel('Latitude', fontsize=15)
gl.right_labels = False
gl.top_labels = False
ax.legend(fontsize=15)
plt.title('N traps: %d; Mean Flux: %0.1f mgC/m^2/d; Mean depth: %dm' % (sum(sel_in_domain),np.mean(POC_in_domain),np.mean(Depth_trap_in_domain)), fontsize=15)
plt.savefig('../Plots/an24/POC_trap_locations_an24.pdf', dpi=200)
plt.close()


########################################################################################################################
#Here I extract only 600m traps
########################################################################################################################
Depth_max_trap = 600
sel_in_domain=(lon>=lon0)&(lon<=lon1)&(lat>=lat0)&(lat<=lat1)&(Depth_trap<=Depth_max_trap)
POC_in_domain=POC[sel_in_domain]




########################################################################################################################
########################################################################################################################
#I load the flux data calculated with our BGC Argo float
########################################################################################################################
########################################################################################################################

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
day0=datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime(2021,10,18)        # final date for the carbon budget calculation
depthf=600                      # final depth
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
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
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
Flux_filtered_depthf=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_filtered_depthf = np.append(Flux_filtered_depthf, np.mean(z[sel_depthf]) )

########################################################################################################################
#I plot flux from our data and from literature
########################################################################################################################
fs=10
width, height = 0.78, 0.75
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.15, width, height])
plt.boxplot([Flux_filtered_depthf,POC_in_domain])
plt.ylabel('POC Flux (mgC/m$^2/d$)', fontsize=fs)
plt.title('Argo Flux calculated between\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day), fontsize=9)
plt.xticks([1,2],['BGC Argo\n6903095 ','Mouw et al.'], fontsize=fs)
plt.savefig('../Plots/an24/01_POCFlux_Argo_vs_Mouw_%d%02d%02dto%d%02d%02d_an24.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()
