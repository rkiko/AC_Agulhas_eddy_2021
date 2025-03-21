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
chla_float_integrated=data_an12['chla_float_integrated']
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

#######################################################################
# I calculate the difference (in terms of percentage) between the mean chl inside
# and outside the eddy (chl which was measured by satellite).
#######################################################################
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum
sel_tmp= np.array(Date_Num_Eddy)>=day0_insideEddy
tmp = np.mean( chl_inside_mean[sel_tmp]-chl_outside_mean[sel_tmp] )
tmp = np.mean( chl_inside_mean[sel_tmp]/chl_outside_mean[sel_tmp] )

#######################################################################
# I plot
#######################################################################

width, height = 0.8, 0.5
set_ylim_lower, set_ylim_upper = -1,1.2#anom_meanE_meanOut.min()*1.1,chl_inside_mean.max()*1.1
fig = plt.figure(1, figsize=(13,4))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_Eddy.min(), Date_Num_Eddy.max()+10))
plt.plot(Date_Num_Eddy,chl_inside_mean,'r',label='Mean Eddy Chl')
plt.plot(Date_Num_Eddy,chl_outside_mean,'g',label='Mean Chl outside')
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
plt.legend(ncol=2,fontsize=12)
plt.ylabel('Eddy Chlorophyll (mg/m$^3$)', fontsize=15)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an20/EddyChl_vs_time_an20.pdf' ,dpi=200)
plt.close()


day_start_eddy_merging=np.array([2021,8,1])
day_end_eddy_merging=np.array([2021,8,11])
day_start_eddy_merging=matlab_datenum(day_start_eddy_merging)-matlab_datenum(1950,1,1)
day_end_eddy_merging=matlab_datenum(day_end_eddy_merging)-matlab_datenum(1950,1,1)
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
# I exclude the 46 element cos it is the profile outside the eddy
sel = [True for i in range(chla_float_integrated.size)]
sel[46]=False
y=chla_float_integrated[sel]
x=Date_Num_float[sel]

fig = plt.figure(1, figsize=(7,3.5))
ax = fig.add_axes([0.12, 0.4, width, height], ylim=(0, y.max()*1.1), xlim=(x.min(), x.max()-27))
plt.plot(x,y,'green')
plt.ylim(ax.get_ylim()[0],ax.get_ylim()[1])
plt.vlines(day_start_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
plt.vlines(day_end_eddy_merging, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color='k')
# I set xticks
nxticks = 10
xticks = np.linspace(x.min(), x.max()-27, nxticks)
xticklabels = []
for i in xticks:
    date_time_obj = date_reference + datetime.timedelta(days=i)
    xticklabels.append(date_time_obj.strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=10)
plt.ylabel('Integrated\nchlorophyll a (mg/m$^2$)', fontsize=10)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an20/IntegratedChl_vs_time_an20.pdf' ,dpi=200)
plt.close()


#######################################################################
# I save the integrated chl concentration values for the latex document
#######################################################################
# from write_latex_data import write_latex_data
# from matlab_datevec import matlab_datevec
# from matlab_datenum import matlab_datenum
# filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
# argument = 'integrated_chl_202104'
# arg_value=np.round(np.mean(y[0]))
# write_latex_data(filename,argument,'%d' % arg_value)
#
# i=47;print(matlab_datevec(x[i]+matlab_datenum(1950,1,1)),y[i])
#
# argument = 'integrated_chl_20210821'
# arg_value=np.round(np.mean(y[47]))
# write_latex_data(filename,argument,'%d' % arg_value)









