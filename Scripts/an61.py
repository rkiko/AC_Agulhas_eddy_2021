import os
import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
import random
from pathlib import Path
home = str(Path.home())
#globals().clear()
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

filename='6903095_Sprof_all.nc'

#######################################################################
# I load the data
#######################################################################
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
Date_Num_float = Date_Num + matlab_datenum(1950,1,1)
pres = np.array(ds.variables['PRES'])
psal = np.array(ds.variables['PSAL'])
temp = np.array(ds.variables['TEMP'])

#######################################################################
# I plot TEMP vs. PSAL for each profile
#######################################################################

# I load sel_insideEddy: for each profile, it tells me whether I am inside the eddy or not
sel_insideEddy = [False for i in range(lon.size)]
a_file = open("%s/an16/data_an16.pkl" % storedir, "rb")
data_an16 = pickle.load(a_file)
tmp=data_an16['sel_insideEddy']
a_file.close()
sel_insideEddy[0:tmp.size]=tmp.copy()

width, height = 0.65, 0.65

fig = plt.figure(1, figsize=(3, 3))
ax = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))

i=0
for i in range(0, pres.shape[0]):
    x = psal[i,:]
    y = temp[i,:]
    sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
    x=x[sel];y=y[sel]
    Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
    plot1=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])

cbar = plt.colorbar(plot1)
nxticks = 10
xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
cbar.set_ticks(xticks)
cbar.set_ticklabels(xticklabels)

plt.xlabel('Salinity (PSU)', fontsize=10)
plt.ylabel('Temperature (°C)', fontsize=10)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an61/TempVsPsal_an61.pdf' , dpi=200)
plt.close()

#Profiles in random order
a=np.array(range(0, pres.shape[0]))
random.shuffle(a)

fig = plt.figure(1, figsize=(3, 3))
ax = fig.add_axes([0.25, 0.2, width, height], ylim=(0, 23), xlim=(34, 36))
i=0
for i in a:
    x = psal[i,:]
    y = temp[i,:]
    sel = (~np.isnan(x)) & (x!=99999) & (~np.isnan(y)) & (y!=99999)
    x=x[sel];y=y[sel]
    Date_Num_tmp = np.tile(Date_Num_float[i], x.size)
    plot1=plt.scatter(x,y,c=Date_Num_tmp,s=0.01,vmin=Date_Num_float.min(),vmax=Date_Num_float[59])

cbar = plt.colorbar(plot1)
nxticks = 10
xticks = np.linspace(Date_Num_float.min(), Date_Num_float[59], nxticks)
xticklabels = []
for i in xticks:
    xticklabels.append('%02d-%02d' % (matlab_datevec(i)[2], matlab_datevec(i)[1]))
cbar.set_ticks(xticks)
cbar.set_ticklabels(xticklabels)

plt.xlabel('Salinity (PSU)', fontsize=10)
plt.ylabel('Temperature (°C)', fontsize=10)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an61/TempVsPsal_random_an61.pdf' , dpi=200)
plt.close()


