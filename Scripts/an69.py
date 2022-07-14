import os,sys
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pandas as pd
import pickle
from CallCppy import Chl_download
from CallCppy import ellipsoid
from scipy.interpolate import griddata
from matplotlib import path
from pathlib import Path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum

#######################################################################################################################
############### Loading Data
#######################################################################################################################


#######################################################################
# I load the eddy center and contours
#######################################################################
filename_eddy1Data='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv' % home
filename_xVMax_eddy1='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy1_lon_an64m.csv' % home
filename_yVMax_eddy1='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy1_lat_an64m.csv' % home
data_eddy1=pd.read_csv(filename_eddy1Data, sep=',', header=0)
data_xVMax1=pd.read_csv(filename_xVMax_eddy1, sep=',', header=0)
data_yVMax1=pd.read_csv(filename_yVMax_eddy1, sep=',', header=0)
Date_Num_Eddy1=data_eddy1['Datenum']
lonEddy1=data_eddy1['Lon']
latEddy1=data_eddy1['Lat']
radius_Vmax1=data_eddy1['Rad_max']   # mean radius_Vmax  = 60.25 km
radius_Out1=data_eddy1['Rad_out']    # mean radius_Out  = 84.60 km
lonVmax1=data_xVMax1.values[:,:]
latVmax1=data_yVMax1.values[:,:]

filename_eddy2Data='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy2.csv' % home
filename_xVMax_eddy2='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy2_lon_an64m.csv' % home
filename_yVMax_eddy2='%s/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy2_lat_an64m.csv' % home
data_eddy2=pd.read_csv(filename_eddy2Data, sep=',', header=0)
data_xVMax2=pd.read_csv(filename_xVMax_eddy2, sep=',', header=0)
data_yVMax2=pd.read_csv(filename_yVMax_eddy2, sep=',', header=0)
Date_Num_Eddy2=data_eddy2['Datenum']
lonEddy2=data_eddy2['Lon']
latEddy2=data_eddy2['Lat']
radius_Vmax2=data_eddy2['Rad_max']   # mean radius_Vmax  = 53.49 km
radius_Out2=data_eddy2['Rad_out']    # mean radius_Out  = 92.01 km
lonVmax2=data_xVMax2.values[:,:]
latVmax2=data_yVMax2.values[:,:]

radius_frame=np.max( [np.mean(radius_Out1), np.mean(radius_Out2)] )*4/111 #radius of the region around the eddy center (in degrees)


#######################################################################
# I save the radius_frame values for the latex document
#######################################################################
# from write_latex_data import write_latex_data
# filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
# argument = 'radius_frame'
# arg_value=radius_frame
# write_latex_data(filename,argument,'%d' % arg_value)

#######################################################################################################################
############### Calculating Chlorophyll using satellite data
#######################################################################################################################

delta_chl=0.04 #Resolution of the chlorophyll product, which is used also to interpolate it

# I load the date range for the chl dt product
a_file = open("%s/GIT/CallCppy/src/CallCppy/Nrt_dt_date.pkl" % (home), "rb")
data_nrt = pickle.load(a_file)
end_date_chl = data_nrt['end_date_chl']
a_file.close()

#### I initialise the chl arrays
chl_inside_mean1=np.squeeze(np.zeros((1,Date_Num_Eddy1.size)))
chl_inside_max1=np.squeeze(np.zeros((1,Date_Num_Eddy1.size)))
chl_outside_mean1=np.squeeze(np.zeros((1,Date_Num_Eddy1.size)))
chl_inside_and_outside_mean1=np.squeeze(np.zeros((1,Date_Num_Eddy1.size)))
chl_inside_mean2=np.squeeze(np.zeros((1,Date_Num_Eddy2.size)))
chl_inside_max2=np.squeeze(np.zeros((1,Date_Num_Eddy2.size)))
chl_outside_mean2=np.squeeze(np.zeros((1,Date_Num_Eddy2.size)))
chl_inside_and_outside_mean2=np.squeeze(np.zeros((1,Date_Num_Eddy2.size)))

#### I start the loop on the two eddies
i_eddy=0
for i_eddy in range(0,2):
    if i_eddy==0:
        Date_Num_Eddy=Date_Num_Eddy1.copy()
        lonEddy=lonEddy1.copy()
        latEddy=latEddy1.copy()
        lonVmax=lonVmax1.copy()
        latVmax=latVmax1.copy()
    if i_eddy==1:
        Date_Num_Eddy=Date_Num_Eddy2.copy()
        lonEddy=lonEddy2.copy()
        latEddy=latEddy2.copy()
        lonVmax = lonVmax2.copy()
        latVmax = latVmax2.copy()

    #### I start the loop on every day for which I have an eddy contour
    i=73
    for i in range(0,lonVmax.shape[1]):#Date_Num_Eddy1.size):
        #### I define the frame around the eddy center, for which I want to download the satellite chl data
        radius_frame_lon=radius_frame / np.cos(lonEddy[i]*np.pi/180)
        lon0_Chl_Down=lonEddy[i]-radius_frame_lon-0.1
        lon1_Chl_Down=lonEddy[i]+radius_frame_lon+0.1
        lat0_Chl_Down=latEddy[i]-radius_frame-0.1
        lat1_Chl_Down=latEddy[i]+radius_frame+0.1
        day0 =  matlab_datevec(Date_Num_Eddy[i]).astype(int)
        print('Starting download Chl for eddy %d on %d-%d-%d..' % (i_eddy+1,day0[2],day0[1],day0[0]))
        #### I download and load the chlorophyll data
        chl_filename,chl_name,lonname,latname=Chl_download(day0, lonmin=lon0_Chl_Down, lonmax=lon1_Chl_Down, latmin=lat0_Chl_Down, latmax=lat1_Chl_Down)

        print('Finished download Chl for eddy %d on %d-%d-%d' % (i_eddy+1,day0[2],day0[1],day0[0]))
        ds = nc.Dataset(chl_filename)
        lon_chl = np.squeeze(np.array(ds.variables[lonname]))
        lat_chl = np.squeeze(np.array(ds.variables[latname]))
        chl_tmp = np.squeeze(np.array(ds.variables[chl_name]))

        # If the Chl_download is nrt, then it downloads the chlorophyll at global level. Therefore, here I select only the part of the chlorophyll matrix which I need
        if np.floor(Date_Num_Eddy[i])>matlab_datenum(end_date_chl):
            a=lon_chl-lon0_Chl_Down
            idx_lon0=int(np.where(a==a[a<0].max())[0])
            a=lon_chl-lon1_Chl_Down
            idx_lon1=int(np.where(a==a[a>0].min())[0])
            a=lat_chl-lat0_Chl_Down
            idx_lat0=int(np.where(a==a[a<0].max())[0])
            a=lat_chl-lat1_Chl_Down
            idx_lat1=int(np.where(a==a[a>0].min())[0])
            lon_chl=lon_chl[idx_lon0:idx_lon1+1] #The longitude goes from -180 to +180
            lat_chl=lat_chl[idx_lat1:idx_lat0+1] #The latitude goes from +90 to -90
            chl_tmp=chl_tmp[idx_lat1:idx_lat0+1,idx_lon0:idx_lon1+1]

        #### I select the eddy maximal velocity contour for the i-th day (excluding the nan values)
        lonVmaxtmp = lonVmax[:,i]
        latVmaxtmp = latVmax[:,i]
        sel=(~np.isnan(lonVmaxtmp)) & (~np.isnan(latVmaxtmp))
        lonVmaxtmp = lonVmaxtmp[sel]
        latVmaxtmp = latVmaxtmp[sel]
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
        print('Finished interpolation for eddy %d on %d-%d-%d' % (i_eddy+1,day0[2],day0[1],day0[0]))
        #### I exclude the nan values (which, in chl_tmp, appear as -999 values)
        sel = chl_interp>=0
        chl_interp = chl_interp[sel]
        sel_outside_contour_Vmax=sel_outside_contour_Vmax[sel]
        #### I extract the chlorophyll values
        if i_eddy == 0:
            chl_inside_mean1[i] = np.nanmean(chl_interp[~sel_outside_contour_Vmax])
            chl_inside_max1[i] = np.nanmax(chl_interp[~sel_outside_contour_Vmax])
            chl_outside_mean1[i] = np.nanmean(chl_interp[sel_outside_contour_Vmax])
            chl_inside_and_outside_mean1[i] = np.nanmean(chl_interp)
        if i_eddy == 1:
            chl_inside_mean2[i] = np.nanmean(chl_interp[~sel_outside_contour_Vmax])
            chl_inside_max2[i] = np.nanmax(chl_interp[~sel_outside_contour_Vmax])
            chl_outside_mean2[i] = np.nanmean(chl_interp[sel_outside_contour_Vmax])
            chl_inside_and_outside_mean2[i] = np.nanmean(chl_interp)

        #### I remove the chlorophyll file
        os.system('rm %s' % chl_filename)


#######################################################################
# Saving data
#######################################################################

dictionary_data = {"chl_inside_mean1": chl_inside_mean1,"chl_inside_max1": chl_inside_max1, "chl_outside_mean1": chl_outside_mean1, "chl_inside_and_outside_mean1": chl_inside_and_outside_mean1,
                   "chl_inside_mean2": chl_inside_mean2,"chl_inside_max2": chl_inside_max2, "chl_outside_mean2": chl_outside_mean2, "chl_inside_and_outside_mean2": chl_inside_and_outside_mean2,
                   "Date_Num_Eddy1": Date_Num_Eddy1, "lonEddy1": lonEddy1, "latEddy1": latEddy1,
                   "Date_Num_Eddy2": Date_Num_Eddy2, "lonEddy2": lonEddy2, "latEddy2": latEddy2}
a_file = open("%s/an69/data_an69.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()





