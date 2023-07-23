
import numpy as np
import pandas as pd
import os,sys
import netCDF4 as nc
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename_coriolis='6903095_Sprof_all.nc'
########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
import seawater as sw
import gsw
from lin_fit import lin_fit

filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
data=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data['RAWfilename']

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
dens=np.array(data['Potential density [kg/m3]'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

Date_Num = np.r_[0:MaP_POC.size]
for i in Date_Num:
    date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    # datetime.datetime.utcfromtimestamp(Date_Num[i])

# I select the data only in the prescribed period
list_dates = np.sort(np.unique(Date_Num))
list_dates = list_dates[(list_dates >= day0_float) & (list_dates <= dayf_float)]
n_profiles = list_dates.size

diffdepth=np.diff(depth)
print(np.mean(diffdepth[diffdepth>0]))

############### I load Coriolis data with the oxygen information
filename='6903095_Sprof_all.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
temp=np.array(ds.variables['TEMP_ADJUSTED'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])
doxy=np.array(ds.variables['DOXY_ADJUSTED'])

i=0;DateVec=np.zeros((Date_Num.size,6)).astype(int)
for i in range(0,Date_Num.size):
    DateVec[i,:]=matlab_datevec(Date_Num[i]+matlab_datenum(1950,1,1))

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
    temp_qc = np.array(ds.variables['TEMP_QC'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
    pres_qc = np.array(ds.variables['PRES_QC'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
    psal_qc = np.array(ds.variables['PSAL_QC'])
if np.sum(doxy==99999)==doxy.size:
    print('Taking non adjusted oxygen')
    doxy = np.array(ds.variables['DOXY'])
    doxy_qc = np.array(ds.variables['DOXY_QC'])

#I tranform the pressure to depth
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth=np.ones(temp.shape)*99999
depth[mask_depth]=depth_tmp

diffdepth=np.diff(depth_tmp)
print(np.mean(diffdepth[diffdepth>0]))

#I calculate the distance between two consecutive profiles
lon0=lon[0:-1]*np.pi/180.0;lon1=lon[1:]*np.pi/180.0
lat0=lat[0:-1]*np.pi/180.0;lat1=lat[1:]*np.pi/180.0
dist = np.arccos(np.sin(lat0 ) * np.sin(lat1 ) + np.cos(lat0 ) * np.cos(lat1 ) * np.cos(lon1 - lon0 )) * 180 / np.pi
dist_km = dist * 111

print(np.mean(dist_km[0:60]))
print(np.max(dist_km[0:60]))
print(np.min(dist_km[0:60]))