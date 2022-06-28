import os
import netCDF4 as nc
import numpy as np
from CallCppy import stadiumShape
from scipy.interpolate import griddata
import pickle
import datetime
from suntime import Sun
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec


dayf=(2021,9,28) # Last day for which I calculate the solar surface irradiance
#######################################################################
# I load the BGC Argo data
#######################################################################

storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename='6903095_Sprof_old.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))
lon_BGC=np.array(ds.variables['LONGITUDE'])
lat_BGC=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)

#######################################################################
# I load and calculate the solar surface irradiance in the days of the BGC Argo profiles and in a given neighborhood
#######################################################################
filename1_ssi = 'S-OSI_-FRA_-MSG_-DLISSID_____-'
filename2_ssi = 'Z.nc'
ray_neighborhood=0.1
delta0 = 0.01

ssi_float=np.array(([]))
ssi_datenum=np.array(([]))
hour_daylight=np.array(([]))
i=0
for i in range(0,Date_Num.size):
    if Date_Num[i]>matlab_datenum(dayf):   continue
    daytmp = matlab_datevec(Date_Num[i])
    filename_ssi = '%s%d%02d%02d1200%s' % (filename1_ssi,daytmp[0],daytmp[1],daytmp[2],filename2_ssi)
    ds = nc.Dataset("../Data/an43/%s" % (filename_ssi))
    lon_ssi = np.array(ds.variables['lon'])
    lat_ssi = np.array(ds.variables['lat'])
    ssi = np.array(ds.variables['ssi']).copy()
    (lon_ssi_g,lat_ssi_g) = np.meshgrid(lon_ssi,lat_ssi)
    (lons, lats) = stadiumShape(np.array([lon_BGC[i], lon_BGC[i]]), np.array([lat_BGC[i], lat_BGC[i]]), ray_neighborhood, delta0)
    ssi_neighborhood = griddata((lon_ssi_g.reshape(lon_ssi_g.size),lat_ssi_g.reshape(lat_ssi_g.size)), ssi.reshape(ssi.size), (lons, lats))
    ssi_float = np.append(ssi_float,np.nanmean(ssi_neighborhood))
    ssi_datenum = np.append(ssi_datenum,Date_Num[i])
    ds.close()
    sun = Sun(lon_BGC[i], lat_BGC[i])
    abd_sr = sun.get_local_sunrise_time(datetime.date(int(daytmp[0]), int(daytmp[1]), int(daytmp[2])))
    abd_ss = sun.get_local_sunset_time(datetime.date(int(daytmp[0]), int(daytmp[1]), int(daytmp[2])))
    hour_daylight = np.append(hour_daylight , (matlab_datenum(abd_ss)-matlab_datenum(abd_sr))*24 )

dictionary_data = {"ssi_per_hour_float": ssi_float, "ssi_datenum": ssi_datenum-matlab_datenum(1950,1,1),
                   "ssi_matlab_datenum" : ssi_datenum, "hour_daylight" : hour_daylight}
a_file = open("%s/an44/data_an44.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()
