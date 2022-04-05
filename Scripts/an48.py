import os
import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#######################################################################
# Parameters
#######################################################################
dist_thr=5 # How many degrees shall a profile be distant to exclude it
time_lag=5 # How many days shall a profile be distant to exclude it
# For the plot:
width, height = 0.65, 0.65
date_reference = datetime.datetime.strptime("1/1/1950", "%d/%m/%Y")
plotdir = '%s/GIT/AC_Agulhas_eddy_2021/Plots/an48/%02ddegree%02ddays' % (home,dist_thr,time_lag)
if os.path.isdir(plotdir) is False:
    os.system('mkdir %s' % plotdir)
    os.system('mkdir %s/PSALvsPRES' % plotdir)
    os.system('mkdir %s/TEMPvsPRES' % plotdir)
    os.system('mkdir %s/TEMPvsPSAL' % plotdir)

#######################################################################
# I load the BGC Argo data of the reference float (690..)
#######################################################################
filename='6903095_Sprof.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon0=np.array(ds.variables['LONGITUDE'])
lat0=np.array(ds.variables['LATITUDE'])
Date_Num0=np.array(ds.variables['JULD'])
pres0 = np.array(ds.variables['PRES'])
psal0 = np.array(ds.variables['PSAL'])
temp0 = np.array(ds.variables['TEMP'])

#######################################################################
# I list the .nc files
#######################################################################
listfiles = os.listdir('%s/an48/' % (storedir))

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
# I start the loop on each profile of the Argo float
#######################################################################

i=0
for i in range(0, pres0.shape[0]):

    lontmp=lon0[i]
    lattmp=lat0[i]
    daytmp=Date_Num0[i]

    #TEMPvsPRES
    x = temp0[i, :]
    y = pres0[i, :]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    fig = plt.figure(1, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 1000), xlim=(0, 23))
    plt.scatter(x, y, c='b', s=0.1,zorder=10)
    plt.plot(x, y, 'b', linewidth=0.25,zorder=11)
    plt.xlabel('Temperature (°C)', fontsize=10)
    plt.ylabel('Pressure (dbar)', fontsize=10)
    plt.gca().invert_yaxis()
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)

    #PSALvsPRES
    x = psal0[i,:]
    y = pres0[i,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    fig = plt.figure(2, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 1000), xlim=(34, 36))
    plt.scatter(x, y, c='b', s=0.1,zorder=10)
    plt.plot(x, y, 'b', linewidth=0.25,zorder=11)
    plt.xlabel('Salinity (PSU)', fontsize=10)
    plt.ylabel('Pressure (dbar)', fontsize=10)
    plt.gca().invert_yaxis()
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)

    #TEMPvsPSAL
    x = psal0[i,:]
    y = temp0[i,:]
    sel = (~np.isnan(x)) & (x != 99999) & (~np.isnan(y)) & (y != 99999)
    x = x[sel];y = y[sel]
    fig = plt.figure(3, figsize=(3, 3))
    ax = fig.add_axes([0.25, 0.15, width, height], ylim=(0, 23), xlim=(34, 36))
    plt.scatter(x, y, c='b', s=0.1,zorder=10)
    plt.plot(x, y, 'b', linewidth=0.25,zorder=11)
    plt.xlabel('Salinity (PSU)', fontsize=10)
    plt.ylabel('Temperature (°C)', fontsize=10)
    plt.gca().invert_yaxis()
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)


    #######################################################################
    # I load the data
    #######################################################################
    n_profiles_ext = 0
    j=0
    for j in range(0,listfiles.__len__()):
        filename=listfiles[j]
        if filename=='.DS_Store':   continue
        ds = nc.Dataset('%s/an48/%s' % (storedir,filename))

        lon=np.array(ds.variables['LONGITUDE'])
        lat=np.array(ds.variables['LATITUDE'])
        Date_Num=np.array(ds.variables['TIME'])

        # if 'PSAL' not in ds.variables:
        #     os.system('rm %s/an48/%s' % (storedir,filename))
        #     print('removed %s' % filename)
        #     continue
        #
        # if 'TEMP' not in ds.variables:
        #     os.system('rm %s/an48/%s' % (storedir,filename))
        #     print('removed %s' % filename)
        #     continue

        if 'PRES_ADJUSTED' not in ds.variables:
            pres = np.array(ds.variables['PRES'])
        else:
            pres = np.array(ds.variables['PRES_ADJUSTED'])

        if 'PSAL_ADJUSTED' not in ds.variables:
            psal = np.array(ds.variables['PSAL'])
        else:
            psal = np.array(ds.variables['PSAL_ADJUSTED'])

        if 'TEMP_ADJUSTED' not in ds.variables:
            temp = np.array(ds.variables['TEMP'])
        else:
            temp = np.array(ds.variables['TEMP_ADJUSTED'])

        sel = np.where( (abs(lon-lontmp)<dist_thr) & (abs(lat-lattmp)<dist_thr) & (abs(Date_Num-daytmp)<time_lag) )[0]
        if sel.size > 0:
            n_profiles_ext = n_profiles_ext + sel.size
            k=0
            for k in range(0,sel.size):
                x = temp[sel[k], :]
                y = pres[sel[k], :]
                z = psal[sel[k], :]
                sel2 = (~np.isnan(x)) & (x != 99999) & (x > -2*10**-9) & (~np.isnan(y)) & (y != 99999) & (y > -2*10**-9) & (~np.isnan(z)) & (z != 99999) & (z > -2*10**-9)
                x = x[sel2];y = y[sel2];z = z[sel2]
                plt.figure(1)
                plt.scatter(x, y, c='gray', s=0.1)
                plt.plot(x, y, 'gray', linewidth=0.25)
                plt.figure(2)
                plt.scatter(z, y, c='gray', s=0.1)
                plt.plot(z, y, 'gray', linewidth=0.25)
                plt.figure(3)
                plt.scatter(z, x, c='gray', s=0.1)
                plt.plot(z, x, 'gray', linewidth=0.25)

    date_time_obj = date_reference + datetime.timedelta(days=Date_Num0[i])
    inoutEddy = 'inside'
    if sel_insideEddy[i] is False:  inoutEddy = 'outside'

    plt.figure(1)
    plt.title('Profile %d, %02d-%02d-%d, %02d:%02d\nFloat %s eddy, %d ext. profiles\n less than %dkm, within %ddays' % (
    int(i + 1), date_time_obj.day, date_time_obj.month, date_time_obj.year, date_time_obj.hour, date_time_obj.minute,
    inoutEddy,n_profiles_ext,dist_thr*111,time_lag), fontsize=10)
    plt.savefig('%s/TEMPvsPRES/TempVsPres_Profile%03d_%02ddegree%02ddays_an48.pdf' % (plotdir,int(i+1),dist_thr,time_lag), dpi=200)
    plt.close()

    plt.figure(2)
    plt.title('Profile %d, %02d-%02d-%d, %02d:%02d\nFloat %s eddy, %d ext. profiles\n less than %dkm, within %ddays' % (
    int(i + 1), date_time_obj.day, date_time_obj.month, date_time_obj.year, date_time_obj.hour, date_time_obj.minute,
    inoutEddy,n_profiles_ext,dist_thr*111,time_lag), fontsize=10)
    plt.savefig('%s/PSALvsPRES/PsalVsPres_Profile%03d_%02ddegree%02ddays_an48.pdf' % (plotdir,int(i+1),dist_thr,time_lag), dpi=200)
    plt.close()

    plt.figure(3)
    plt.title('Profile %d, %02d-%02d-%d, %02d:%02d\nFloat %s eddy, %d ext. profiles\n less than %dkm, within %ddays' % (
    int(i + 1), date_time_obj.day, date_time_obj.month, date_time_obj.year, date_time_obj.hour, date_time_obj.minute,
    inoutEddy,n_profiles_ext,dist_thr*111,time_lag), fontsize=10)
    plt.savefig('%s/TEMPvsPSAL/TempVsPsal_Profile%03d_%02ddegree%02ddays_an48.pdf' % (plotdir,int(i+1),dist_thr,time_lag), dpi=200)
    plt.close()




