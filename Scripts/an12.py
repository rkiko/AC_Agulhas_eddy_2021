import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4 as nc
import seawater as sw
import gsw
import pandas as pd
import pickle
from paruvpy import mixed_layer_depth
from CallCppy import Chl_download
from CallCppy import ellipsoid
from pathlib import Path
from scipy.interpolate import griddata
from matplotlib import path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#######################################################################################################################
############### Loading Data
#######################################################################################################################

#######################################################################
# I load the Coriolis data
#######################################################################

# To update the data, please run:
# os.system("python Download_BGC_variables.py")
filename='6903095_Sprof.nc'

ds = nc.Dataset('%s/%s' % (storedir,filename))
lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])

Date_Num=np.array(ds.variables['JULD'])
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
Date_Vec=np.zeros([Date_Num.size,6])
for i in range(0,Date_Num.size):
    date_time_obj = date_reference + datetime.timedelta(days=Date_Num[i])
    Date_Vec[i,0]=date_time_obj.year;Date_Vec[i,1]=date_time_obj.month;Date_Vec[i,2]=date_time_obj.day
    Date_Vec[i,3]=date_time_obj.hour;Date_Vec[i,4]=date_time_obj.minute;Date_Vec[i,5]=date_time_obj.second

Date_Vec=Date_Vec.astype(int)

#I load the variables
temp=np.array(ds.variables['TEMP_ADJUSTED'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])
chla=np.array(ds.variables['CHLA_ADJUSTED'])

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
if np.sum(chla==99999)==chla.size:
    print('Taking non adjusted chlorophyll-a')
    chla = np.array(ds.variables['CHLA'])

#I tranform the pressure to depth
mask_depth=pres!=99999 #I select only valid values
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_depth]
pres_tmp=pres[mask_depth]
depth_tmp=sw.eos80.dpth(pres_tmp, lat_tmp)
depth=np.ones(temp.shape)*99999
depth[mask_depth]=depth_tmp

#I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
#salinity and temperature first
mask_dens=np.logical_and(pres!=99999,temp!=99999,psal!=99999) # I exclude the points with value = 99999
lat_tmp=np.tile(lat,[pres.shape[1],1]).T
lon_tmp=np.tile(lon,[pres.shape[1],1]).T
lat_tmp=lat_tmp[mask_dens]
lon_tmp=lon_tmp[mask_dens]
pres_tmp=pres[mask_dens]
psal_tmp=psal[mask_dens]
temp_tmp=temp[mask_dens]
abs_psal_tmp=gsw.SA_from_SP(psal_tmp, pres_tmp, lon_tmp, lat_tmp) # I compute absolute salinity
cons_tmp=gsw.CT_from_t(abs_psal_tmp, temp_tmp, pres_tmp)          # I compute conservative temperature
temp_tmp=gsw.density.sigma0(abs_psal_tmp, cons_tmp)
dens=np.ones(temp.shape)*99999
dens[mask_dens]=temp_tmp+1000

#######################################################################
# I load the eddy center and contours
#######################################################################
filename_eddyData='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_data_TOEddies.csv' % home
filename_xVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Xcon_max.csv' % home
filename_yVMax='%s/GIT/AC_Agulhas_eddy_2021/Data/an12/BE_cyclone_TOEddies_Ycon_max.csv' % home
data_eddy=pd.read_csv(filename_eddyData, sep=',', header=0)
data_xVMax=pd.read_csv(filename_xVMax, sep=',', header=None)
data_yVMax=pd.read_csv(filename_yVMax, sep=',', header=None)
lonEddy=data_eddy['x_centroid']
latEddy=data_eddy['y_centroid']
radius_Vmax=data_eddy['Rmax']   # mean radius_Vmax  = 60 km
radius_Out=data_eddy['Rout']    # mean radius_Out  = 85 km
radius_frame=np.mean(radius_Out)*4/111 #radius of the region around the eddy center (in degrees)
lonVmax=data_xVMax.values[:,:]
latVmax=data_yVMax.values[:,:]
Date_Num_Eddy=data_eddy['date_num']
DateTime_Eddy=[]
i=Date_Num_Eddy[0]
for i in Date_Num_Eddy:
    DateTime_Eddy=np.append(DateTime_Eddy,datetime.datetime.fromordinal(int(i)-366))

#######################################################################################################################
############### Calculating Chlorophyll using BGC Argo data
#######################################################################################################################

#######################################################################
# Calculating Chlorophyll along the BGC Argo float trajectory
#######################################################################

chla_float_mean=np.squeeze(np.zeros((1,chla.shape[0])))
chla_float_max=np.squeeze(np.zeros((1,chla.shape[0])))

i=0
for i in range(0,chla.shape[0]):
    chla_tmp=chla[i,:]
    depth_tmp=depth[i,:]
    temp_tmp=temp[i,:]
    # I exclude nan values
    sel_non_nan=(chla_tmp!=99999)&(depth_tmp!=99999)&(temp_tmp!=99999)
    chla_tmp=chla_tmp[sel_non_nan];temp_tmp=temp_tmp[sel_non_nan];depth_tmp=depth_tmp[sel_non_nan]
    mld,_ = mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
    # I select only the chl values inside the mixed layer
    sel_in_ML=depth_tmp<=mld
    chla_tmp=chla_tmp[sel_in_ML]
    chla_float_mean[i]=np.mean(chla_tmp)
    chla_float_max[i]=np.max(chla_tmp)

#######################################################################################################################
############### Calculating Chlorophyll using satellite data
#######################################################################################################################

delta_chl=0.04 #Resolution of the chlorophyll product, which is used also to inteprolate it
#### I initialise the chl arrays
chl_inside_mean=np.squeeze(np.zeros((1,Date_Num_Eddy.size)))
chl_inside_max=np.squeeze(np.zeros((1,Date_Num_Eddy.size)))
chl_outside_mean=np.squeeze(np.zeros((1,Date_Num_Eddy.size)))
chl_inside_and_outside_mean=np.squeeze(np.zeros((1,Date_Num_Eddy.size)))

#### I start the loop on every day for which I have an eddy contour
i=0
for i in range(0,Date_Num_Eddy.size):
    #### I define the frame around the eddy center, for which I want to download the satellite chl data
    radius_frame_lon=radius_frame / np.cos(lonEddy[i]*np.pi/180)
    lon0_Chl_Down=lonEddy[i]-radius_frame_lon-0.1
    lon1_Chl_Down=lonEddy[i]+radius_frame_lon+0.1
    lat0_Chl_Down=latEddy[i]-radius_frame-0.1
    lat1_Chl_Down=latEddy[i]+radius_frame+0.1
    day0 =  np.array((DateTime_Eddy[i].year,DateTime_Eddy[i].month,DateTime_Eddy[i].day))
    #### I download and load the chlorophyll data
    chl_filename,chl_name,lonname,latname=Chl_download(day0, lonmin=lon0_Chl_Down, lonmax=lon1_Chl_Down, latmin=lat0_Chl_Down, latmax=lat1_Chl_Down)
    ds = nc.Dataset(chl_filename)
    lon_chl = np.squeeze(np.array(ds.variables[lonname]))
    lat_chl = np.squeeze(np.array(ds.variables[latname]))
    chl_tmp = np.squeeze(np.array(ds.variables[chl_name]))
    #### I select the eddy maximal velocity contour for the i-th day (excluding the nan values)
    lonVmaxtmp = lonVmax[:,i]
    latVmaxtmp = latVmax[:,i]
    sel=(~np.isnan(lonVmaxtmp)) & (~np.isnan(latVmaxtmp))
    lonVmaxtmp = lonVmaxtmp[sel]
    latVmaxtmp = latVmaxtmp[sel]
    #plt.plot(lonVmaxtmp,latVmaxtmp,'r.');plt.plot(lonVmaxtmp,latVmaxtmp,'r')
    # I create a circular cloud frame around the eddy center
    lons, lats = ellipsoid(lonEddy[i], latEddy[i], radius_frame_lon, radius_frame, delta_chl)
    #plt.plot(lons,lats,'bo');plt.plot(lonVmaxtmp,latVmaxtmp,'r.');plt.plot(lonVmaxtmp,latVmaxtmp,'r')
    #### I select only the points of the circular cloud frame which are not inside the eddy maximal velocity contour
    contour_Vmax=np.concatenate((lonVmaxtmp.reshape(lonVmaxtmp.size,1),latVmaxtmp.reshape(latVmaxtmp.size,1)),axis=1)
    p = path.Path(contour_Vmax)
    frame_points = np.concatenate((lons.reshape(lons.size,1),lats.reshape(lats.size,1)),axis=1)
    sel_outside_contour_Vmax = ~p.contains_points(frame_points)
    #plt.plot(lons,lats,'bo');plt.plot(lonVmaxtmp,latVmaxtmp,'r.');plt.plot(lonVmaxtmp,latVmaxtmp,'r');plt.plot(lons[sel_outside_contour_Vmax],lats[sel_outside_contour_Vmax],'ko');
    #### I interpolate the chlorophyll values for each point of the circular cloud frame
    lon_chl_g, lat_chl_g = np.meshgrid(lon_chl, lat_chl)
    lon_chl_g, lat_chl_g, chl_tmp = np.squeeze(lon_chl_g.reshape(lon_chl_g.size, 1)), np.squeeze(lat_chl_g.reshape(lat_chl_g.size, 1)), np.squeeze(chl_tmp.reshape(chl_tmp.size, 1))
    chl_interp = griddata((lon_chl_g, lat_chl_g), chl_tmp, (lons, lats))
    #### I exclude the nan values (which, in chl_tmp, appear as -999 values)
    sel = chl_interp>=0
    chl_interp = chl_interp[sel]
    sel_outside_contour_Vmax=sel_outside_contour_Vmax[sel]
    #### I extract the chlorophyll values
    chl_inside_mean[i] = np.nanmean(chl_interp[~sel_outside_contour_Vmax])
    chl_inside_max[i] = np.nanmax(chl_interp[~sel_outside_contour_Vmax])
    chl_outside_mean[i] = np.nanmean(chl_interp[sel_outside_contour_Vmax])
    chl_inside_and_outside_mean[i] = np.nanmean(chl_interp)
    #### I remove the chlorophyll file
    os.system('rm %s' % chl_filename)


#######################################################################
# Saving data
#######################################################################

dictionary_data = {"chla_float_mean": chla_float_mean, "chla_float_max": chla_float_max, "chl_inside_mean": chl_inside_mean,
                   "chl_inside_max": chl_inside_max, "chl_outside_mean": chl_outside_mean, "chl_inside_and_outside_mean": chl_inside_and_outside_mean,
                   "Date_Num_Eddy": Date_Num_Eddy, "lonEddy": lonEddy, "latEddy": latEddy,
                   "lon_float": lon, "lat_float": lat, "Date_Num_float": Date_Num, "Date_Vec_float": Date_Vec}
a_file = open("%s/an12/data_an12.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()





