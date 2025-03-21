import os
import ftplib
import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4 as nc
import seawater as sw
import gsw
from pathlib import Path
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
home = str(Path.home())
#globals().clear()
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#usr="anonymous"
#pwd="alberto.baudena@gmail.com"
#HOSTNAME="usgodae.org"
#ftp_server = ftplib.FTP(HOSTNAME, usr,pwd)
#ftp_server.cwd('pub/outgoing/argo/dac/coriolis/6903095')
#ftp_server.dir()
#ftp_server.dir('profiles')
#ftp_server.cwd('profiles')

filename='6903095_BRtraj.nc'
filename='6903095_Rtraj.nc'
filename='BR6903095_001.nc'
filename='6903095_Sprof_old.nc'

#os.chdir(storedir)
#ftp_server.retrbinary("RETR " + filename, open('%s' % filename, 'wb').write)
#os.chdir(actualdir)
#ftp_server.quit()

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

#Standard variables
temp=np.array(ds.variables['TEMP_ADJUSTED'])
temp_qc=np.array(ds.variables['TEMP_ADJUSTED_QC'])
temp_qc_profile=np.array(ds.variables['PROFILE_TEMP_QC'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
pres_qc=np.array(ds.variables['PRES_ADJUSTED_QC'])
pres_qc_profile=np.array(ds.variables['PROFILE_PRES_QC'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])
psal_qc=np.array(ds.variables['PSAL_ADJUSTED_QC'])
psal_qc_profile=np.array(ds.variables['PROFILE_PSAL_QC'])

#BGC Variables
chla=np.array(ds.variables['CHLA_ADJUSTED'])
chla_qc=np.array(ds.variables['CHLA_ADJUSTED_QC'])
chla_qc_profile=np.array(ds.variables['PROFILE_CHLA_QC'])
doxy=np.array(ds.variables['DOXY_ADJUSTED'])
doxy_qc=np.array(ds.variables['DOXY_ADJUSTED_QC'])
doxy_qc_profile=np.array(ds.variables['PROFILE_DOXY_QC'])
bbp700=np.array(ds.variables['BBP700_ADJUSTED'])
bbp700_qc=np.array(ds.variables['BBP700_ADJUSTED_QC'])
bbp700_qc_profile=np.array(ds.variables['PROFILE_BBP700_QC'])

#If adjusted values are not available yet, I take the non adjusted ones
if np.sum(temp==99999)==temp.size:
    print('Taking non adjusted temperature')
    temp = np.array(ds.variables['TEMP'])
    temp_qc = np.array(ds.variables['TEMP_QC'])
if np.sum(pres==99999)==pres.size:
    print('Taking non adjusted pressure')
    pres = np.array(ds.variables['PRES'])
    pres_qc = np.array(ds.variables['PRES_QC'])
if np.sum(psal==99999)==psal.size:
    print('Taking non adjusted salinity')
    psal = np.array(ds.variables['PSAL'])
    psal_qc = np.array(ds.variables['PSAL_QC'])
if np.sum(chla==99999)==chla.size:
    print('Taking non adjusted chlorophyll-a')
    chla = np.array(ds.variables['CHLA'])
    chla_qc = np.array(ds.variables['CHLA_QC'])
if np.sum(doxy==99999)==doxy.size:
    print('Taking non adjusted oxygen')
    doxy = np.array(ds.variables['DOXY'])
    doxy_qc = np.array(ds.variables['DOXY_QC'])
if np.sum(bbp700==99999)==bbp700.size:
    print('Taking non adjusted bbp700')
    bbp700 = np.array(ds.variables['BBP700'])
    bbp700_qc = np.array(ds.variables['BBP700_QC'])

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
dens_tmp=gsw.density.sigma0(abs_psal_tmp, cons_tmp)
dens=np.ones(temp.shape)*99999
dens[mask_dens]=dens_tmp+1000

#Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
#isopycnal, given by reference_isopycnal
reference_isopycnal=1027
parameter=doxy
Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+117])
doxy_isopycnal=np.array([]);Date_Num_isopycnal=np.array([])
i=0
for i in range(0,doxy.shape[0]):
    sel = (parameter[i, :] != 99999) & (dens[i, :] != 99999)
    z=doxy[i,sel];y=dens[i,sel]
    for iy in range(0,y.size-1):
        if np.logical_and(Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1], y[iy] <= reference_isopycnal < y[iy + 1]):
            dist = (reference_isopycnal - y[iy]) / (y[iy + 1] - y[iy])
            doxy_tmp= z[iy] + (z[iy + 1] - z[iy]) * dist
            doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
            Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])

#Here I interpolate, so that the slope gives me the respiration rate (in micromol/kg per day)
from AC_Agulhas_eddy_2021.Scripts.lin_fit import lin_fit
(interpol,_,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)

#Here I plot
width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(12, 8))
ax = fig.add_axes([0.12, 0.2, width, height])
plot3 = plt.scatter(Date_Num_isopycnal, doxy_isopycnal, c='black')
plot4 = plt.plot(np.linspace(Date_Num_limit[0], Date_Num_limit[1], 20),
                 np.linspace(Date_Num_limit[0] * interpol.slope + interpol.intercept,
                             Date_Num_limit[1] * interpol.slope + interpol.intercept, 20), c='black')
plt.text(Date_Num_limit[0]+(Date_Num_limit[1]-Date_Num_limit[0])*0.1,doxy_isopycnal.min()+(doxy_isopycnal.max()-doxy_isopycnal.min())*0.2,'R$^2$=%.2f, p=%0.2e' % (interpol.rvalue**2,interpol.pvalue), fontsize=16)
plt.ylabel('Dissolved oxygen ($\mu$mol/kg)', fontsize=18)
plt.title('Doxy time series along isopycnal %.1f kg/m$^3$, slope: %0.3f $\mu$mol/kg/d fit signif: %s' % (reference_isopycnal,interpol.slope,signif_label), fontsize=15)
# I set xticks
nxticks = 10
xticks = np.linspace(Date_Num_limit[0], Date_Num_limit[1], nxticks)
xticklabels = []
for i in xticks:
    date_time_obj = date_reference + datetime.timedelta(days=i)
    xticklabels.append(date_time_obj.strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90, fontsize=12)
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an06/OxygenTimeSeries_isopycnal%0.1f_an06.pdf' % reference_isopycnal, dpi=200)






