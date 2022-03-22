import os
import numpy as np
import netCDF4 as nc
import numpy as np
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datenum import matlab_datenum
from matlab_datevec import matlab_datevec

day0=(2021,4,13)
dayf=(2021,9,28)
lonmin,lonmax = 5,18
latmin,latmax = -40, -30
ndays = int(matlab_datenum(dayf) - matlab_datenum(day0)+1)
filename1 = 'S-OSI_-FRA_-MSG_-DLISSID_____-'
filename2= 'Z.nc'


i=0
for i in range(0,ndays):
    daytmp = matlab_datevec(matlab_datenum (day0)+i )
    filename = '%s%d%02d%02d1200%s' % (filename1,daytmp[0],daytmp[1],daytmp[2],filename2)
    filenameB = '%s%d%02d%02d1200%sB.nc' % (filename1,daytmp[0],daytmp[1],daytmp[2],filename2[:-3])
    ds = nc.Dataset("../Data/an43/%s" % (filename))
    lon2 = np.array(ds.variables['lon'])
    lat2 = np.array(ds.variables['lat'])
    ssib = np.array(ds.variables['ssi'])
    landmaskb = np.array(ds.variables['landmask'])

    # import matplotlib.pyplot as plt
    # (lon2g, lat2g) = np.meshgrid(lon2,lat2)
    # fig = plt.figure(2, figsize=(12, 8))
    # plot1 = plt.contourf(lon2,lat2,ssib)
    # plot2 = plt.plot(lon2g[landmaskb==1],lat2g[landmaskb==1],'k.')
    # cbar = plt.colorbar(plot1)
    # plt.xlim(lonmin, lonmax)
    # plt.ylim(latmin,latmax)

    ilonmin=np.where(np.abs(lon2-lonmin) == (np.abs(lon2-lonmin)).min())[0][0]
    ilonmax=np.where(np.abs(lon2-lonmax) == (np.abs(lon2-lonmax)).min())[0][0]
    ilatmin=np.where(np.abs(lat2-latmin) == (np.abs(lat2-latmin)).min())[0][0]
    ilatmax=np.where(np.abs(lat2-latmax) == (np.abs(lat2-latmax)).min())[0][0]
    lon2 = lon2[ilonmin:ilonmax]
    lat2 = lat2[ilatmin:ilatmax]
    ssib = ssib[ilatmin:ilatmax,ilonmin:ilonmax]
    landmaskb = landmaskb[ilatmin:ilatmax,ilonmin:ilonmax]
    ds.close()
    ds = nc.Dataset("../Data/an43/%s" % (filenameB),mode='w')
    ds.createDimension('lat', lat2.size)
    ds.createDimension('lon', lon2.size)
    lat = ds.createVariable('lat', np.float32, ('lat',))
    lon = ds.createVariable('lon', np.float32, ('lon',))
    ssi = ds.createVariable('ssi', np.float32, ('lat','lon'))
    landmask = ds.createVariable('landmask', np.float32, ('lat','lon'))
    lon[:] = lon2
    lat[:] = lat2
    ssi[:,:] = ssib
    landmaskb[:,:] = landmaskb
    ds.close()
    os.system("rm -f ../Data/an43/%s" % (filename))
    os.system("mv ../Data/an43/%s ../Data/an43/%s" % (filenameB,filename))
