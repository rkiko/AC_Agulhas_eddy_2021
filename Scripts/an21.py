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
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

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

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Date_Time.size]
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
radius_Out=data_eddy['Rout']
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
Date_Num_Eddy_4float=np.array(Date_Num_Eddy[sel_Satel2Float[sel_eddyCont_missing_days]])

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
# I save the radius Vmax values for the latex document
#######################################################################
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from write_latex_data import write_latex_data
from matlab_datevec import matlab_datevec

Date_Vec_Eddy=np.zeros((Date_Num_Eddy.size,6))
for i in range(0,Date_Num_Eddy.size):
    Date_Vec_Eddy[i,:] = matlab_datevec(Date_Num_Eddy[i])

Date_Vec_Eddy = Date_Vec_Eddy.astype(int)
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
argument = 'radius_Vmax_20201018'
arg_value=radius_Vmax[0]
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'radius_Vmax_20210413'
arg_value=radius_Vmax[175]
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'max_radius_Vmax'
arg_value=radius_Vmax.max()
write_latex_data(filename,argument,'%d' % arg_value)
arg_value=Date_Vec_Eddy[np.where(radius_Vmax == radius_Vmax.max()),:]
write_latex_data(filename,'max_radius_Vmax_date','%d August' % arg_value[0][0][2])

#######################################################################
# I plot
#######################################################################

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = 0,150
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_Eddy_4float[0], Date_Num_Eddy_4float.max()))
plt.plot(Date_Num_Eddy,radius_Vmax,'b',label='Mean Chl')
plt.plot(Date_Num_Eddy,radius_Out,'r--',label='Chl Anomaly')
plt.scatter(Date_Num_Eddy_4float,dist_km,c='k',label='Chl Anomaly',marker='*')
# I set xticks
nxticks = 10
xticks = np.linspace(Date_Num_Eddy_4float[0], Date_Num_Eddy.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.fromordinal(int(i)-366).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
# plt.legend(fontsize=12)
plt.ylabel('Radius (km)', fontsize=15)
ax.text(-0.075, 1.05, 'a', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an21/EddyRadiusAndDist_vs_time_an21.pdf' ,dpi=200)
plt.close()





