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


dayf=(2021,9,28) # Last day for which I calculate the euphotic depth
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
# I load and calculate the euphotic layer depth in the days of the BGC Argo profiles and in a given neighborhood
#######################################################################
year=2021
filename1_zeu='A'
filename2_zeu = 'L3m_8D_ZLEE_Zeu_lee_9km.nc'
ray_neighborhood=0.1
delta0 = 0.01
list_index_day0=np.r_[1:305:8]
list_index_dayf=np.r_[8:305:8]

zeu_float=np.array(([]))
zeu_datenum=np.array(([]))
i=0
for i in range(0,Date_Num.size):
    if Date_Num[i]>matlab_datenum(dayf):   continue
    index_day = int(Date_Num[i]-matlab_datenum(2020,12,31))
    index_day = index_day - list_index_day0
    index_day[index_day < 0] = 99999
    index_day = np.where(index_day == index_day.min())[0][0]
    index0=list_index_day0[index_day]
    indexf=list_index_dayf[index_day]
    filename_zeu = '%s%d%03d%d%03d.%s' % (filename1_zeu, year, index0, year, indexf, filename2_zeu)
    ds = nc.Dataset("../Data/an35/8D/%s" % (filename_zeu))
    lon_zeu = np.array(ds.variables['lon'])
    lat_zeu = np.array(ds.variables['lat'])
    zeu = np.array(ds.variables['Zeu_lee']).copy()
    zeu[zeu<-10]=np.nan
    (lon_zeu_g,lat_zeu_g) = np.meshgrid(lon_zeu,lat_zeu)
    (lons, lats) = stadiumShape(np.array([lon_BGC[i], lon_BGC[i]]), np.array([lat_BGC[i], lat_BGC[i]]), ray_neighborhood, delta0)
    zeu_neighborhood = griddata((lon_zeu_g.reshape(lon_zeu_g.size),lat_zeu_g.reshape(lat_zeu_g.size)), zeu.reshape(zeu.size), (lons, lats))
    zeu_float = np.append(zeu_float,np.nanmean(zeu_neighborhood))
    zeu_datenum = np.append(zeu_datenum,Date_Num[i])
    ds.close()

dictionary_data = {"zeu_float": zeu_float, "zeu_datenum": zeu_datenum-matlab_datenum(1950,1,1),
                   "zeu_matlab_datenum" : zeu_datenum}
a_file = open("%s/an36/data_an36.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()
