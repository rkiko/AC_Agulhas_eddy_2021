import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

########################################################################################################################
#Time period that I plot
########################################################################################################################
day0=datetime(2021,4,13)        # starting date
dayf=datetime(2021,9,24)        # final date
day0_float = calendar.timegm(day0.timetuple())
dayf_float = calendar.timegm(dayf.timetuple())
ndays = (dayf - day0).days  # number of days
delta_bin=20 #thickness (in meters) of the bin I used to interpolate the psd slope

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
diffPSD_slop=np.array(data['diffPSD slope'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:depth.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]

########################################################################################################################
#I filter the data
########################################################################################################################

diffPSD_slop_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=diffPSD_slop[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        iy=0
        for iy in range(0,y2.size):
            sel=(y2>=y2[iy]-delta_bin*0.5)&(y2<y2[iy]+delta_bin*0.5)
            z[iy]=np.mean(z[sel])
        diffPSD_slop_filtered = np.concatenate((diffPSD_slop_filtered, z))
        Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
        depth_filtered = np.concatenate((depth_filtered, y2))

########################################################################################################################
#I interpolate the data
########################################################################################################################
# I define the x and y arrays for the MiP+MaP+bbp interpolation
x_filtered = np.linspace(Date_Num_filtered.min(), Date_Num_filtered.max(), int(np.round(ndays/4)))
y_filtered = np.linspace(depth_filtered.min(), depth_filtered.max(), 100)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
diffPSD_slop_interp = griddata((Date_Num_filtered, depth_filtered), diffPSD_slop_filtered,(x_filtered_g, y_filtered_g), method="nearest")

########################################################################################################################
#I plot the filtered data
########################################################################################################################

width, height = 0.8, 0.68
set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
fig = plt.figure(1, figsize=(7,3.5))
ax = fig.add_axes([0.12, 0.25, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
parameter_plot=diffPSD_slop_interp.copy()
# parameter_plot[parameter_plot<0]=0
# parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
plt.gca().invert_yaxis()
# I draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('diff PSD slope', fontsize=10)
plt.ylabel('Depth (m)', fontsize=10)
plt.title('Filtered diff PSD slope, delta bin: %d m' % delta_bin, fontsize=10)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=7)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an26/01Filtered_diffPSD_slop_bin%dm_an26.pdf' % delta_bin,dpi=200)
plt.close()








########################################################################################################################
#I take the unfiltered data
########################################################################################################################

diffPSD_slop_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=diffPSD_slop[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        # z = savgol_filter(z, 5, 1)
        iy=0
        for iy in range(0,y2.size):
            sel=(y2>=y2[iy]-delta_bin*0.5)&(y2<y2[iy]+delta_bin*0.5)
            z[iy]=np.mean(z[sel])
        diffPSD_slop_filtered = np.concatenate((diffPSD_slop_filtered, z))
        Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
        depth_filtered = np.concatenate((depth_filtered, y2))

########################################################################################################################
#I interpolate the data
########################################################################################################################
# I define the x and y arrays for the MiP+MaP+bbp interpolation
x_filtered = np.linspace(Date_Num_filtered.min(), Date_Num_filtered.max(), int(np.round(ndays/4)))
y_filtered = np.linspace(depth_filtered.min(), depth_filtered.max(), 100)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
diffPSD_slop_interp = griddata((Date_Num_filtered, depth_filtered), diffPSD_slop_filtered,(x_filtered_g, y_filtered_g), method="nearest")

########################################################################################################################
#I plot the unfiltered data
########################################################################################################################

width, height = 0.8, 0.68
set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
fig = plt.figure(1, figsize=(7,3.5))
ax = fig.add_axes([0.12, 0.25, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_filtered.min(), Date_Num_filtered.max()))
parameter_plot=diffPSD_slop_interp.copy()
# parameter_plot[parameter_plot<0]=0
# parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
plt.gca().invert_yaxis()
# I draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('diff PSD slope', fontsize=10)
plt.ylabel('Depth (m)', fontsize=10)
plt.title('Non Filtered diff PSD slope, delta bin: %d m' % delta_bin, fontsize=10)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=7)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an26/00NonFiltered_diffPSD_slop_bin%dm_an26.pdf' % delta_bin,dpi=200)
plt.close()
