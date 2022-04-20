import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from datetime import datetime
import netCDF4 as nc
import gsw
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
plotdir = '%s/GIT/AC_Agulhas_eddy_2021/Plots/an52/' % (home)
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

########################################################################
#Parameters
########################################################################
day_release0=matlab_datenum(2021,4,13)
day_release1=matlab_datenum(2021,4,15)

########################################################################
# I load the CTD data
########################################################################
filename_alpha='%s/GIT/AC_Agulhas_eddy_2021/Data/an52/SO283_CTD.tab' % home

data=pd.read_csv(filename_alpha, sep='\t', header=110)
lon=data.Longitude
lat=data.Latitude
event=data['Event']
depth=data['Depth water [m] (CTD, SEA-BIRD SBE 9 plus)']
Date_Time=data['Date/Time']
Press=data['Press [dbar] (CTD, SEA-BIRD SBE 9 plus)']
Temp1=data['Temp [°C] (CTD, SEA-BIRD SBE 9 plus)']
Temp2=data['Temp [°C] (2, CTD, SEA-BIRD SBE 9 plus)']
Sal1=data['Sal (PSU, CTD, SEA-BIRD SBE 9 plus)']
Sal2=data['Sal (2, PSU, CTD, SEA-BIRD SBE 9 plus)']
Doxy1=data['DO [ml/l] (CTD, SEA-BIRD SBE 9 plus)']
Doxy2=data['DO [ml/l] (2, CTD, SEA-BIRD SBE 9 plus)']

Date_Num=np.squeeze( np.zeros((Date_Time.size,1)) )
event_list = np.r_[0:Date_Time.size]
i=0
for i in range(0,Date_Time.size):
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = matlab_datenum(date_time_obj.year,date_time_obj.month,date_time_obj.day,date_time_obj.hour,date_time_obj.minute,date_time_obj.second)
    event_list[i]  = event[i].split('_')[1].split('-')[0]

list_dates=np.sort(np.unique(Date_Num))
event_list_unique= np.sort(np.unique(event_list))

########################################################################
# I load the BGC Argo float data
########################################################################

filename='6903095_Sprof_all.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon0=np.array(ds.variables['LONGITUDE'])
lat0=np.array(ds.variables['LATITUDE'])
Date_Num0=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
pres0 = np.array(ds.variables['PRES'])
psal0 = np.array(ds.variables['PSAL'])
temp0 = np.array(ds.variables['TEMP'])
doxy0 = np.array(ds.variables['DOXY_ADJUSTED']) #micromol/kg

########################################################################
# I convert the oxygen of the BGC Argo float from micromol/kg to ml/l
########################################################################
mask_doxy=doxy0 > 99990
liter_per_mole0=22.414 #approximation at 0 degree
liter_per_mole=liter_per_mole0*(temp0+273.15)/273.15 #approximation at 0 degree
doxy0 = doxy0*liter_per_mole #microl/kg
doxy0 = doxy0/1000 #ml/kg

# To transform from ml/kg to ml/l, I compute the potential density
mask_dens = np.logical_and(pres0 != 99999, temp0 != 99999, psal0 != 99999)  # I exclude the points with value = 99999
lat0_tmp = np.tile(lat0, [pres0.shape[1], 1]).T
lon0_tmp = np.tile(lon0, [pres0.shape[1], 1]).T
lat0_tmp = lat0_tmp[mask_dens]
lon0_tmp = lon0_tmp[mask_dens]
pres0_tmp = pres0[mask_dens]
psal0_tmp = psal0[mask_dens]
temp0_tmp = temp0[mask_dens]
abs_psal0_tmp = gsw.SA_from_SP(psal0_tmp, pres0_tmp, lon0_tmp, lat0_tmp)  # I compute absolute salinity
cons_tmp = gsw.CT_from_t(abs_psal0_tmp, temp0_tmp, pres0_tmp)  # I compute conservative temperature
dens0_tmp = gsw.density.sigma0(abs_psal0_tmp, cons_tmp)
dens0 = np.ones(temp0.shape) * 99999
dens0[mask_dens] = dens0_tmp + 1000

# I transform from ml/kg to ml/l
doxy0 = doxy0 * dens0/1000
doxy0[mask_doxy|~mask_dens] = 99999

########################################################################
# I select the events which were taken time_lag day close to the float release event
########################################################################

sel = (Date_Num>day_release0)&(Date_Num<day_release1)
event_close_ctd = event_list[sel]
event_close_ctd_unique= np.sort(np.unique(event_close_ctd))

########################################################################
# I plot
########################################################################

width, height = 0.7, 0.75
fig = plt.figure(1, figsize=(3, 3))
ax1 = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 500))#, xlim=(0, 23))
fig = plt.figure(2, figsize=(3, 3))
ax2 = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 500))#, xlim=(0, 23))
fig = plt.figure(3, figsize=(3, 3))
ax3 = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 500))#, xlim=(0, 10)
fig = plt.figure(4, figsize=(3, 3))
ax4 = fig.add_axes([0.25, 0.15, width, height])#, ylim=(0, 500))#, xlim=(0, 23))
for j in range(0,6):
    daytmp = matlab_datevec(Date_Num0[j])

    plt.figure(1)
    x = temp0[j, :]; y = pres0[j,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x,y,label='%d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2],daytmp[3],daytmp[4],abs(lat0[j]),lon0[j]) )

    plt.figure(2)
    x = psal0[j, :]; y = pres0[j,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x,y,label='%d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2],daytmp[3],daytmp[4],abs(lat0[j]),lon0[j]) )

    plt.figure(3)
    x = doxy0[j, :]; y = pres0[j,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x,y,label='%d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2],daytmp[3],daytmp[4],abs(lat0[j]),lon0[j]) )

    plt.figure(4)
    x = psal0[j, :]; y = temp0[j,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x,y,label='%d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2],daytmp[3],daytmp[4],abs(lat0[j]),lon0[j]) )

i = event_close_ctd_unique[0]
for i in event_close_ctd_unique:
    sel = event_list == i

    lon_tmp = np.array(lon[sel])
    lat_tmp = np.array(lat[sel])
    depth_tmp = np.array(depth[sel])
    Date_Num_tmp = np.array(Date_Num[sel])
    daytmp = matlab_datevec(Date_Num_tmp[0])
    Press_tmp = np.array(Press[sel])
    Temp1_tmp = np.array(Temp1[sel])
    Temp2_tmp = np.array(Temp2[sel])
    Sal1_tmp = np.array(Sal1[sel])
    Sal2_tmp = np.array(Sal2[sel])
    Doxy1_tmp = np.array(Doxy1[sel])
    Doxy2_tmp = np.array(Doxy2[sel])

    plt.figure(1)
    x = Temp1_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD1, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )
    x = Temp2_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD2, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )

    plt.figure(2)
    x = Sal1_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD1, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )
    x = Sal2_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD2, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )

    plt.figure(3)
    x = Doxy1_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD1, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )
    x = Doxy2_tmp;y = Press_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD2, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )

    plt.figure(4)
    x = Sal1_tmp;y = Temp1_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD1, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )
    x = Sal2_tmp;y = Temp2_tmp
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    plt.plot(x, y, label='CTD2, %d Apr, %02dh%02d\n%0.2fS,%0.2fE' % (daytmp[2], daytmp[3], daytmp[4],abs(lat_tmp[0]),lon_tmp[0]) )

plt.figure(1)
plt.xlabel('Temperature (°C)', fontsize=10)
plt.ylabel('Pressure (dbar)', fontsize=10)
plt.gca().invert_yaxis()
plt.legend(fontsize=5)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
ax1.text(-0.215, 1.05, 'a', transform=ax1.transAxes,fontsize=14, fontweight='bold', va='top', ha='right')
plt.savefig('%s/TempVsPres_an52.pdf' % (plotdir), dpi=200)
plt.close()

plt.figure(2)
plt.xlabel('Salinity (Psu)', fontsize=10)
plt.ylabel('Pressure (dbar)', fontsize=10)
plt.gca().invert_yaxis()
plt.legend(fontsize=5)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
ax2.text(-0.215, 1.05, 'b', transform=ax2.transAxes,fontsize=14, fontweight='bold', va='top', ha='right')
plt.savefig('%s/SalVsPres_an52.pdf' % (plotdir), dpi=200)
plt.close()

plt.figure(3)
plt.xlabel('Oxygen (ml/l)', fontsize=10)
plt.ylabel('Pressure (dbar)', fontsize=10)
plt.gca().invert_yaxis()
plt.legend(fontsize=5)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
ax3.text(-0.215, 1.05, 'd', transform=ax3.transAxes,fontsize=14, fontweight='bold', va='top', ha='right')
plt.savefig('%s/DoxyVsPres_an52.pdf' % (plotdir), dpi=200)
plt.close()

plt.figure(4)
plt.xlabel('Salinity (Psu)', fontsize=10)
plt.ylabel('Temperature (°C)', fontsize=10)
plt.legend(fontsize=5)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
ax4.text(-0.215, 1.05, 'c', transform=ax4.transAxes,fontsize=14, fontweight='bold', va='top', ha='right')
plt.savefig('%s/TempVsSal_an52.pdf' % (plotdir), dpi=200)
plt.close()


