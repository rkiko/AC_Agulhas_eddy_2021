import os
import netCDF4 as nc
import numpy as np
from CallCppy import stadiumShape
from scipy.interpolate import griddata
import pickle
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum


dayf=(2021,9,28) # Last day for which I calculate the attenuation coefficient k490
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
# I load and calculate the attenuation coefficient k490 in the days of the BGC Argo profiles and in a given neighborhood
#######################################################################
year=2021
filename1_Kd_490='A'
filename2_Kd_490 = 'L3m_8D_KD490_Kd_490_9km.nc'
ray_neighborhood=0.1
delta0 = 0.01
list_index_day0=np.r_[1:305:8]
list_index_dayf=np.r_[8:305:8]

Kd_490_float=np.array(([]))
Kd_490_datenum=np.array(([]))
i=0
for i in range(0,Date_Num.size):
    if Date_Num[i]>matlab_datenum(dayf):   continue
    index_day = int(Date_Num[i]-matlab_datenum(2020,12,31))
    index_day = index_day - list_index_day0
    index_day[index_day < 0] = 99999
    index_day = np.where(index_day == index_day.min())[0][0]
    index0=list_index_day0[index_day]
    indexf=list_index_dayf[index_day]
    filename_Kd_490 = '%s%d%03d%d%03d.%s' % (filename1_Kd_490, year, index0, year, indexf, filename2_Kd_490)
    ds = nc.Dataset("../Data/an41/8D/%s" % (filename_Kd_490))
    lon_Kd_490 = np.array(ds.variables['lon'])
    lat_Kd_490 = np.array(ds.variables['lat'])
    Kd_490 = np.array(ds.variables['Kd_490']).copy()
    Kd_490[Kd_490<-10]=np.nan
    (lon_Kd_490_g,lat_Kd_490_g) = np.meshgrid(lon_Kd_490,lat_Kd_490)
    (lons, lats) = stadiumShape(np.array([lon_BGC[i], lon_BGC[i]]), np.array([lat_BGC[i], lat_BGC[i]]), ray_neighborhood, delta0)
    Kd_490_neighborhood = griddata((lon_Kd_490_g.reshape(lon_Kd_490_g.size),lat_Kd_490_g.reshape(lat_Kd_490_g.size)), Kd_490.reshape(Kd_490.size), (lons, lats))
    Kd_490_float = np.append(Kd_490_float,np.nanmean(Kd_490_neighborhood))
    Kd_490_datenum = np.append(Kd_490_datenum,Date_Num[i])
    ds.close()

dictionary_data = {"Kd_490_float": Kd_490_float, "Kd_490_datenum": Kd_490_datenum-matlab_datenum(1950,1,1),
                   "Kd_490_matlab_datenum" : Kd_490_datenum}
a_file = open("%s/an42/data_an42.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()
