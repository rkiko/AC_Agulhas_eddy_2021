import calendar
import numpy as np
import os
import pandas as pd
from datetime import datetime
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

########################################################################################################################
#
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
lon_float=np.array(data['Longitude'][sel_filename])
lat_float=np.array(data['Latitude'][sel_filename])
Date_Time_float=np.array(data['Date_Time'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num_float_calendar=np.r_[0:Date_Time_float.size]
Date_Num_float=np.r_[0:Date_Time_float.size].astype(float)
for i in Date_Num_float_calendar:
    date_time_obj = datetime.strptime(Date_Time_float[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num_float[i] = matlab_datenum(date_time_obj.year,date_time_obj.month,date_time_obj.day,date_time_obj.hour,date_time_obj.minute,date_time_obj.second)
    Date_Num_float_calendar[i] = calendar.timegm(date_time_obj.timetuple())

########################################################################################################################
# I load the eddy radius and position
########################################################################################################################
filename_eddyData='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies.csv' % home
data_eddy=pd.read_csv(filename_eddyData, sep=',', header=0)
lonEddy=data_eddy['x_centroid']
latEddy=data_eddy['y_centroid']
Date_Num_Eddy=data_eddy['date_num']
radius_Vmax=data_eddy['Rmax']
radius_Out=data_eddy['Rout']

########################################################################################################################
# I do a loop on every profile
########################################################################################################################
list_dates=np.sort(np.unique( Date_Num_float))[0:61]
dist_float_eddy_km = np.array([])
radius_Vmax_floatDays = np.array([])
radius_Out_floatDays = np.array([])
i=0
for i in range (0,list_dates.size):
    Date_Num_float_tmp = np.floor( list_dates[i] )
    i_float = np.where( np.floor(Date_Num_float) == Date_Num_float_tmp )[0][0]
    # I take the closest day in which I have an eddy contour for Date_Num_float_tmp
    tmp = abs( Date_Num_Eddy - Date_Num_float_tmp )
    i_eddy = np.where( tmp == tmp.min() )[0][0]

    # I extract the float and eddy parameters for the day Date_Num_float_tmp
    lon_float_tmp = lon_float [i_float]
    lat_float_tmp = lat_float [i_float]
    lon_eddy_tmp = lonEddy [i_eddy]
    lat_eddy_tmp = latEddy [i_eddy]
    radius_Vmax_tmp = radius_Vmax [i_eddy]
    radius_Out_tmp = radius_Out [i_eddy]

    # I calculate the distance
    dist = np.arccos(np.sin(lat_float_tmp * np.pi / 180) * np.sin(lat_eddy_tmp * np.pi / 180) + np.cos(
        lat_float_tmp * np.pi / 180) * np.cos(lat_eddy_tmp * np.pi / 180) * np.cos(
        (lon_eddy_tmp - lon_float_tmp) * np.pi / 180)) * 180 / np.pi
    dist_km = dist * 111

    # I append the data
    dist_float_eddy_km = np.append(dist_float_eddy_km, dist_km)
    radius_Vmax_floatDays = np.append(radius_Vmax_floatDays, radius_Vmax_tmp)
    radius_Out_floatDays = np.append(radius_Out_floatDays, radius_Out_tmp)

########################################################################################################################
# I save
########################################################################################################################
if list_dates.size==69:
    dictionary_data = {"Date_Num_float": list_dates,"dist_float_eddy_km": dist_float_eddy_km,
                       "radius_Vmax_floatDays": radius_Vmax_floatDays,"radius_Out_floatDays": radius_Out_floatDays}
    a_file = open("%s/GIT/AC_Agulhas_eddy_2021/Data/an57/data_an57.pkl" % (home), "wb")
    pickle.dump(dictionary_data, a_file)
    a_file.close()


########################################################################################################################
# I plot
########################################################################################################################
width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = 0,150
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(list_dates[0], list_dates.max()))
plt.plot(Date_Num_Eddy,radius_Vmax,'b')
# plt.plot(list_dates,radius_Vmax_floatDays,'bo')
plt.plot(Date_Num_Eddy,radius_Out,'r--')
# plt.plot(list_dates,radius_Out_floatDays,'ro')
plt.scatter(list_dates,dist_float_eddy_km,c='k',marker='*')
# I set xticks
nxticks = 10
xticks = np.linspace(list_dates[0], list_dates.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.fromordinal(int(i)-366).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.ylabel('Radius (km)', fontsize=15)
ax.text(-0.075, 1.05, 'a', transform=ax.transAxes,fontsize=34, fontweight='bold', va='top', ha='right') # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.xticks(rotation=90, fontsize=14)
plt.savefig('../Plots/an57/Eddy_float_distance_plus radius_an57.pdf' ,dpi=200)
plt.close()

########################################################################################################################
#
########################################################################################################################







