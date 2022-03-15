import os
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

year=2021
day0=103 #2021,4,13
dayf=266 #2021,9,23

download_webiste='https://oceandata.sci.gsfc.nasa.gov/ob/getfile'
filename0='A'
filename2 = 'L3m_DAY_ZLEE_Zeu_lee_9km.nc'
appkey='a69b30815ad5f6246b09769c97b6cd244f751f6d'

index_day=day0
for index_day in range(day0,dayf+1):
    filename = '%s%d%03d.%s' % (filename0,year,index_day,filename2)
    url_filename_download = '%s/%s?appkey=%s' % (download_webiste, filename,appkey)
    os.system("python obdaac_download.py -v %s" % (url_filename_download))
    os.system("mv %s ../Data/an35/" % (filename))

import netCDF4 as nc
import numpy as np
lonmin,lonmax = 5,18
latmin,latmax = -40, -30

index_day=day0
for index_day in range(day0,dayf+1):
    filename = '%s%d%03d.%s' % (filename0,year,index_day,filename2)
    filenameB = '%s%d%03d.%sB.nc' % (filename0,year,index_day,filename2[:-3])
    ds = nc.Dataset("../Data/an35/%s" % (filename))
    lon2 = np.array(ds.variables['lon'])
    lat2 = np.array(ds.variables['lat'])
    zeu2 = np.array(ds.variables['Zeu_lee'])
    ilonmin=np.where(np.abs(lon2-lonmin) == (np.abs(lon2-lonmin)).min())[0][0]
    ilonmax=np.where(np.abs(lon2-lonmax) == (np.abs(lon2-lonmax)).min())[0][0]
    ilatmin=np.where(np.abs(lat2-latmin) == (np.abs(lat2-latmin)).min())[0][0]
    ilatmax=np.where(np.abs(lat2-latmax) == (np.abs(lat2-latmax)).min())[0][0]
    lon2 = lon2[ilonmin:ilonmax]
    lat2 = lat2[ilatmax:ilatmin]
    zeu2 = zeu2[ilatmax:ilatmin,ilonmin:ilonmax]
    ds.close()
    ds = nc.Dataset("../Data/an35/%s" % (filenameB),mode='w')
    ds.createDimension('lat', lat2.size)
    ds.createDimension('lon', lon2.size)
    lat = ds.createVariable('lat', np.float32, ('lat',))
    lon = ds.createVariable('lon', np.float32, ('lon',))
    zeu = ds.createVariable('Zeu_lee', np.float32, ('lat','lon'))
    lon[:] = lon2
    lat[:] = lat2
    zeu[:,:] = zeu2
    ds.close()
    os.system("rm -f ../Data/an35/%s" % (filename))
    os.system("mv ../Data/an35/%s ../Data/an35/%s" % (filenameB,filename))
