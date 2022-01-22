import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import pickle
from datetime import date
from pathlib import Path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#######################################################################################################################
############### Loading Data
#######################################################################################################################
a_file = open("%s/an12/data_an12.pkl" % storedir, "rb")
data_an12 = pickle.load(a_file)
chla_float_mean=data_an12['chla_float_mean']
chla_float_max=data_an12['chla_float_max']
lon_float=data_an12['lon_float']
lat_float=data_an12['lat_float']
Date_Num_float=data_an12['Date_Num_float']
Date_Vec_float=data_an12['Date_Vec_float']
chl_inside_mean=data_an12['chl_inside_mean']
chl_inside_max=data_an12['chl_inside_max']
chl_outside_mean=data_an12['chl_outside_mean']
chl_inside_and_outside_mean=data_an12['chl_inside_and_outside_mean']
Date_Num_Eddy=data_an12['Date_Num_Eddy']
DateTime_Eddy=data_an12['DateTime_Eddy']
lonEddy=data_an12['lonEddy']
latEddy=data_an12['latEddy']
a_file.close()

# ### I load eddy contour data
# filename_xVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Xcon_max.csv' % home
# filename_yVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Ycon_max.csv' % home
# data_xVMax=pd.read_csv(filename_xVMax, sep=',', header=None)
# data_yVMax=pd.read_csv(filename_yVMax, sep=',', header=None)
# lonVmax=data_xVMax.values[:,:]
# latVmax=data_yVMax.values[:,:]
#
# ### I load eddy radius
# filename_eddyData='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies.csv' % home
# data_eddy=pd.read_csv(filename_eddyData, sep=',', header=0)
# radius_Vmax=data_eddy['Rmax']
# del data_eddy

#######################################################################
# I select the data only in the period when the BGC Argo float was inside the eddy
#######################################################################
day0_insideEddy=[2021,4,13]
dayf_insideEddy=[2021,9,27]
day0_insideEddy=date.toordinal(date(day0_insideEddy[0], day0_insideEddy[1], day0_insideEddy[2]))+366
dayf_insideEddy=date.toordinal(date(dayf_insideEddy[0], dayf_insideEddy[1], dayf_insideEddy[2]))+366


###Satellite data: here the difference is that I consider also the days before the BGC argo entered the eddy
### I don't need to filter the eddy contour data as they are updated and contain only eddy contours of the period in which
### BGC Argo was inside the eddy
sel_insideEddy = [False for i in range(lonEddy.size)]
i=0
for i in range(0,lonEddy.size):
    if (Date_Num_Eddy[i]<=dayf_insideEddy):
        sel_insideEddy[i] = True

chl_inside_mean=chl_inside_mean[sel_insideEddy]
chl_inside_max=chl_inside_max[sel_insideEddy]
chl_outside_mean=chl_outside_mean[sel_insideEddy]
chl_inside_and_outside_mean=chl_inside_and_outside_mean[sel_insideEddy]
Date_Num_Eddy=Date_Num_Eddy[sel_insideEddy]
DateTime_Eddy=DateTime_Eddy[sel_insideEddy]
lonEddy=lonEddy[sel_insideEddy]
latEddy=latEddy[sel_insideEddy]
anom_meanE_meanOut=chl_inside_mean-chl_outside_mean




width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = -1,1#anom_meanE_meanOut.min()*1.1,chl_inside_mean.max()*1.1
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_Eddy.min(), Date_Num_Eddy.max()+10))
plt.plot(Date_Num_Eddy,chl_inside_mean,'r',label='Mean Chl')
plt.plot(Date_Num_Eddy,anom_meanE_meanOut,'b',label='Chl Anomaly')
plt.vlines([day0_insideEddy,dayf_insideEddy],set_ylim_lower, set_ylim_upper,colors='black',label='Float inside eddy',linewidth=3,linestyles='dashed')
plt.hlines([set_ylim_lower+0.05,set_ylim_upper-0.05],day0_insideEddy, dayf_insideEddy,colors='black',linewidth=3,linestyles='dashed')
# I set xticks
nxticks = 10
xticks = np.linspace(Date_Num_Eddy.min(), Date_Num_Eddy.max(), nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append(datetime.datetime.fromordinal(int(i)-366).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=14)
plt.legend(fontsize=12)
plt.ylabel('Eddy Chlorophyll (mg/m$^3$)', fontsize=15)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an20/EddyChl_vs_time_an20.pdf' ,dpi=200)
plt.close()












