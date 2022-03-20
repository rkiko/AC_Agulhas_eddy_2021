import os
import numpy as np
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

year=2021
day0=103 #2021,4,13
dayf=271 #2021,9,28

list_index_day0=np.r_[1:305:8]
list_index_dayf=np.r_[8:305:8]
index0=day0-list_index_day0
index0[index0<0]=99999
index0=np.where(index0==index0.min())[0][0]
indexf=dayf-list_index_day0
indexf[indexf<0]=99999
indexf=np.where(indexf==indexf.min())[0][0]

list_index_day0=list_index_day0[index0:indexf+1]
list_index_dayf=list_index_dayf[index0:indexf+1]

download_webiste='https://oceandata.sci.gsfc.nasa.gov/ob/getfile'
filename1='A'
filename2 = 'L3m_8D_KD490_Kd_490_9km.nc'
appkey='a69b30815ad5f6246b09769c97b6cd244f751f6d'

i=0
for i in range(0,list_index_day0.size):
    filename = '%s%d%03d%d%03d.%s' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2)
    url_filename_download = '%s/%s?appkey=%s' % (download_webiste, filename,appkey)
    os.system("python obdaac_download.py -v %s" % (url_filename_download))
    os.system("mv %s ../Data/an41/8D/" % (filename))
    # To make the plot
    # import netCDF4 as nc
    # import numpy as np
    # lonmin, lonmax = 5, 18
    # latmin, latmax = -40, -30
    # filename = 'A20210972021104.L3m_8D_KD490_Kd_490_9km.nc'
    # ds = nc.Dataset("../Data/an41/8D/%s" % (filename))
    # lon2 = np.array(ds.variables['lon'])
    # lat2 = np.array(ds.variables['lat'])
    # Kd_490 = np.array(ds.variables['Kd_490'])
    # Kd_490b=Kd_490.copy()
    # Kd_490b[Kd_490b<-10]=np.nan
    # Kd_490b[Kd_490b>1000]=np.nan
    # Kd_490b[Kd_490b>0.25]=0.25
    # import matplotlib.pyplot as plt
    # fig = plt.figure(1, figsize=(12, 8))
    # plot1 = plt.contourf(lon2,lat2,Kd_490b)
    # plt.xlim(lonmin,lonmax)
    # plt.ylim(latmin,latmax)
    # cbar = plt.colorbar(plot1)
    # ds.close()





import netCDF4 as nc
import numpy as np
lonmin,lonmax = 5,18
latmin,latmax = -40, -30

i=0
for i in range(0,list_index_day0.size):
    filename = '%s%d%03d%d%03d.%s' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2)
    filenameB = '%s%d%03d%d%03d.%sB.nc' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2[:-3])
    ds = nc.Dataset("../Data/an41/8D/%s" % (filename))
    lon2 = np.array(ds.variables['lon'])
    lat2 = np.array(ds.variables['lat'])
    Kd_490b = np.array(ds.variables['Kd_490'])
    ilonmin=np.where(np.abs(lon2-lonmin) == (np.abs(lon2-lonmin)).min())[0][0]
    ilonmax=np.where(np.abs(lon2-lonmax) == (np.abs(lon2-lonmax)).min())[0][0]
    ilatmin=np.where(np.abs(lat2-latmin) == (np.abs(lat2-latmin)).min())[0][0]
    ilatmax=np.where(np.abs(lat2-latmax) == (np.abs(lat2-latmax)).min())[0][0]
    lon2 = lon2[ilonmin:ilonmax]
    lat2 = lat2[ilatmax:ilatmin]
    Kd_490b = Kd_490b[ilatmax:ilatmin,ilonmin:ilonmax]
    ds.close()
    ds = nc.Dataset("../Data/an41/8D/%s" % (filenameB),mode='w')
    ds.createDimension('lat', lat2.size)
    ds.createDimension('lon', lon2.size)
    lat = ds.createVariable('lat', np.float32, ('lat',))
    lon = ds.createVariable('lon', np.float32, ('lon',))
    Kd_490 = ds.createVariable('Kd_490', np.float32, ('lat','lon'))
    lon[:] = lon2
    lat[:] = lat2
    Kd_490[:,:] = Kd_490b
    ds.close()
    os.system("rm -f ../Data/an41/8D/%s" % (filename))
    os.system("mv ../Data/an41/8D/%s ../Data/an41/8D/%s" % (filenameB,filename))
