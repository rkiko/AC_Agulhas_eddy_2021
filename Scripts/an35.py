import os
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
filename2 = 'L3m_8D_ZLEE_Zeu_lee_9km.nc'
appkey='a69b30815ad5f6246b09769c97b6cd244f751f6d'

i=0
for i in range(0,list_index_day0.size):
    filename = '%s%d%03d%d%03d.%s' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2)
    url_filename_download = '%s/%s?appkey=%s' % (download_webiste, filename,appkey)
    os.system("python obdaac_download.py -v %s" % (url_filename_download))
    os.system("mv %s ../Data/an35/8D/" % (filename))
    # To make the plot
    # import netCDF4 as nc
    # import numpy as np
    # lonmin, lonmax = 5, 18
    # latmin, latmax = -40, -30
    #
    # filename = 'A20211052021112.L3m_8D_ZLEE_Zeu_lee_9km.nc'
    # ds = nc.Dataset("%s" % (filename))
    # lon2 = np.array(ds.variables['lon'])
    # lat2 = np.array(ds.variables['lat'])
    # zeu2 = np.array(ds.variables['Zeu_lee'])
    # zeu3=zeu2.copy()
    # zeu3[zeu3<-10]=np.nan
    # zeu3[zeu3>1000]=np.nan
    # zeu3[zeu3>150]=150
    # import matplotlib.pyplot as plt
    # fig = plt.figure(1, figsize=(12, 8))
    # plot1 = plt.contourf(lon2,lat2,zeu3)
    # plt.xlim(lonmin,lonmax)
    # plt.ylim(latmin,latmax)
    # cbar = plt.colorbar(plot1)
    # ds.close()
    # filename = '%s%d%03d.%s' % (filename0, year, index_day, filename2)
    # ds = nc.Dataset("../Data/an35/8D/%s" % (filename))
    # lon_zeu = np.array(ds.variables['lon'])
    # lat_zeu = np.array(ds.variables['lat'])
    # zeu = np.array(ds.variables['Zeu_lee']).copy()
    # zeu[zeu<-10]=np.nan
    # zeu[zeu>1000]=np.nan
    # zeu[zeu>150]=150
    # fig2 = plt.figure(2, figsize=(12, 8))
    # plot2 = plt.contourf(lon_zeu,lat_zeu,zeu)
    # cbar = plt.colorbar(plot2)





import netCDF4 as nc
import numpy as np
lonmin,lonmax = 5,18
latmin,latmax = -40, -30

i=0
for i in range(0,list_index_day0.size):
    filename = '%s%d%03d%d%03d.%s' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2)
    filenameB = '%s%d%03d%d%03d.%sB.nc' % (filename1,year,list_index_day0[i],year,list_index_dayf[i],filename2[:-3])
    ds = nc.Dataset("../Data/an35/8D/%s" % (filename))
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
    ds = nc.Dataset("../Data/an35/8D/%s" % (filenameB),mode='w')
    ds.createDimension('lat', lat2.size)
    ds.createDimension('lon', lon2.size)
    lat = ds.createVariable('lat', np.float32, ('lat',))
    lon = ds.createVariable('lon', np.float32, ('lon',))
    zeu = ds.createVariable('Zeu_lee', np.float32, ('lat','lon'))
    lon[:] = lon2
    lat[:] = lat2
    zeu[:,:] = zeu2
    ds.close()
    os.system("rm -f ../Data/an35/8D/%s" % (filename))
    os.system("mv ../Data/an35/8D/%s ../Data/an35/8D/%s" % (filenameB,filename))
