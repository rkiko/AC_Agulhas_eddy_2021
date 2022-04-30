import numpy as np
import pandas as pd
import os
from scipy.io import loadmat
from scipy.interpolate import griddata
import seawater as sw
import netCDF4 as nc
import datetime
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
home = str(Path.home())
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from CallCppy import ellipsoid
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum

#######################################################################
# Parameters
#######################################################################
list_radius=np.array([0.1,0.25,0.5,1,2]) # The radius of the neighborhood I use to extract the doxy climatology for a given profile location
delta0=0.05

# For the plot:
width, height = 0.65, 0.65
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
plotdir = '%s/GIT/AC_Agulhas_eddy_2021/Plots/an54' % (home)

#######################################################################
#I load the climatological doxy variables
#######################################################################
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
mat = loadmat("../Data/global_DO_field.mat")
doxy_clim = mat['DATAgridded_DO']
lon = mat['xvec'][0]
lat = mat['yvec'][0]
pres = np.squeeze(mat['pvec'])

#######################################################################
# I load the BGC Argo data of the reference float (690..)
#######################################################################
filename='6903095_Sprof_all.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon0=np.array(ds.variables['LONGITUDE'])
lat0=np.array(ds.variables['LATITUDE'])
Date_Num0=np.array(ds.variables['JULD'])
pres0 = np.array(ds.variables['PRES'])
psal0 = np.array(ds.variables['PSAL'])
temp0 = np.array(ds.variables['TEMP'])
doxy0 = np.array(ds.variables['DOXY_ADJUSTED'])

#######################################################################
# I load sel_insideEddy: for each profile, it tells me whether I am inside the eddy or not
#######################################################################
sel_insideEddy = [False for i in range(lon0.size)]
a_file = open("%s/an16/data_an16.pkl" % storedir, "rb")
data_an16 = pickle.load(a_file)
tmp=data_an16['sel_insideEddy']
a_file.close()
sel_insideEddy[0:tmp.size]=tmp.copy()


#######################################################################
#I reduce the climatology to save time
#######################################################################
radius_max=list_radius.max()
lonmin=np.floor(lon0.min()-radius_max)
lonmax=np.ceil(lon0.max()+radius_max)
latmin=np.floor(lat0.min()-radius_max)
latmax=np.ceil(lat0.max()+radius_max)
lonmin = np.where(lon==lonmin)[0][0]
lonmax = np.where(lon==lonmax)[0][0]
latmin = np.where(lat==latmin)[0][0]
latmax = np.where(lat==latmax)[0][0]

lon=lon[lonmin:lonmax+1]
lat=lat[latmin:latmax+1]
doxy_clim = doxy_clim[lonmin:lonmax+1,latmin:latmax+1,:]
lon_g=np.tile(lon,[lat.size,1]).T.reshape((lon.size*lat.size))
lat_g=np.tile(lat,[lon.size,1]).reshape((lon.size*lat.size))

#######################################################################
# I start the loop on each profile of the Argo float
#######################################################################

i=84
for i in range(0, pres0.shape[0]):

    x = doxy0[i, :]
    y = pres0[i, :]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    below1000m_YN = 1
    if y.max() < 1500:    below1000m_YN = 2
    fig = plt.figure(4, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.15, width, height], ylim=(0000, 2000))#, xlim=(180, 230))
    plt.scatter(x, y, c='b', s=0.1,zorder=10)
    plt.plot(x, y, 'b', linewidth=0.25,zorder=11)
    plt.xlabel('Dissolved oxygen ($\mu$mol/kg)', fontsize=10)
    plt.ylabel('Pressure (dbar)', fontsize=10)
    plt.gca().invert_yaxis()
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)

    # I extract the doxy climatology
    lon1=lon0[i]
    lat1=lat0[i]
    j=0
    for j in range(0,list_radius.__len__()):
        radius=list_radius[j]
        lons, lats = ellipsoid(lon1, lat1, radius, radius, delta0)

        x2=np.zeros((pres.size,))
        k=0
        for k in range(0,x2.size):
            doxy_tmp = doxy_clim[:, :, k].reshape(lon_g.size,)
            doxy_interp = griddata((lon_g, lat_g), doxy_tmp, (lons, lats))
            x2[k] = doxy_interp.mean()

        plt.plot(x2, pres, linewidth=0.25, zorder=10-j, label = 'R = %d km' % int(radius*111))

    plt.legend(fontsize=6)
    date_time_obj = date_reference + datetime.timedelta(days=Date_Num0[i])
    inoutEddy = 'inside'
    if sel_insideEddy[i] is False:  inoutEddy = 'outside'
    plt.title('Profile %d, %02d-%02d-%d, %02d:%02d\nFloat %s eddy' % (
        int(i + 1), date_time_obj.day, date_time_obj.month, date_time_obj.year, date_time_obj.hour,
        date_time_obj.minute,inoutEddy), fontsize=10)

    plt.savefig('%s/%02dProfile%02d_an54.pdf' % ( plotdir,below1000m_YN,int(i + 1) ), dpi=200)
    plt.close()
