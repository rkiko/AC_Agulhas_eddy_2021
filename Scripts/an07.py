import calendar
import numpy as np
import os
import pandas as pd
from datetime import timedelta,datetime
from scipy.interpolate import interp1d
import netCDF4 as nc
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory


####################################################################################
#####################ECOPART DATA
####################################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_processed_data.tsv' % home
data=pd.read_csv(filename_ecopart, sep='\t', header=0)
RAWfilename=data.RAWfilename

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data
lon=np.array(data['Longitude'][sel_filename]);lat=np.array(data['Latitude'][sel_filename])
Date_Time=np.array(data['Date_Time'][sel_filename])
RAWfilename=RAWfilename[sel_filename]
pressure_EcoPart=-np.array(data['Pressure [dbar]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:lon.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

####################################################################################
#####################CORIOLIS DATA
####################################################################################

storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename_Coriolis='6903095_Sprof.nc'
ds = nc.Dataset('%s/%s' % (storedir,filename_Coriolis))
lon2=np.array(ds.variables['LONGITUDE'])
lat2=np.array(ds.variables['LATITUDE'])

Date_Num2=np.array(ds.variables['JULD'])
date_reference = datetime.strptime("1/1/1950", "%d/%m/%Y")
Date_Vec=np.zeros([Date_Num.size,6])
for i in range(0,Date_Num2.size):
    date_time_obj = date_reference + timedelta(days=Date_Num2[i])
    Date_Num2[i] = calendar.timegm(date_time_obj.timetuple())
    Date_Vec[i,0]=date_time_obj.year;Date_Vec[i,1]=date_time_obj.month;Date_Vec[i,2]=date_time_obj.day
    Date_Vec[i,3]=date_time_obj.hour;Date_Vec[i,4]=date_time_obj.minute;Date_Vec[i,5]=date_time_obj.second
    # datetime.utcfromtimestamp(Date_Num2[53])

Date_Num2=Date_Num2.astype(int)
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

####################################################################################
#####################LOOP ON ECOPART DATA, PROFILE PER PROFILE
####################################################################################

list_dates=np.sort(np.unique(Date_Num))
if list_dates.size!=pres.shape[0]: raise ValueError("Error: EcoPart and Coriolis data do not have the same number of profiles. Please update them")

i=0;temp_EcoPart=np.squeeze(np.zeros((lon.size,1)));doxy_EcoPart=np.squeeze(np.zeros((lon.size,1)))
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    # I comment all these lines in which I was trying to identify the Coriolis profile from the longitude+latitude of the
    # Ecopart profile since it does not work. This because the longitude and latitude does not match very well.
    # However, I found a new solution (for the moment) which is based selecting the Ecopart profiles progressively according
    # to the date. The Coriolis profiles have also a date which is always incrementing, so they (should) match.
    # RAWfilename_profile=np.unique(RAWfilename[sel])
    # if RAWfilename_profile.size>1: raise ValueError("%d raw_filenames associated to the Ecopart profile carried in data %s" % (RAWfilename_profile.size,datetime.utcfromtimestamp(list_dates[i]).__str__()))
    # lon_profile = np.unique(lon[sel])
    # lat_profile = np.unique(lat[sel])
    # if lon_profile.size>1: raise ValueError("Error: the %s Ecopart profile has %d longitude values" % (RAWfilename_profile[0],lon_profile.size))
    # if lat_profile.size>1: raise ValueError("Error: the %s Ecopart profile has %d latitude values" % (RAWfilename_profile[0],lat_profile.size))
    # idx_lon = np.array(np.where(np.round(lon2*1000)==int(lon_profile*1000)))
    # idx_lat = np.array(np.where(np.round(lat2*1000)==int(lat_profile*1000)))
    # if idx_lon.size>1: raise ValueError("Error: the %s Ecopart profile is associated with %d different Coriolis longitude profiles" % (RAWfilename_profile[0],idx_lon.size))
    # if idx_lat.size>1: raise ValueError("Error: the %s Ecopart profile is associated with %d different Coriolis latitude profiles" % (RAWfilename_profile[0],idx_lat.size))
    # if idx_lon!=idx_lat: raise ValueError("Error: the %s Ecopart profile associated with the Coriolis profile using longitude is different from the Ecopart profile associated with the Coriolis profile using latitude (%d and %d)" % (RAWfilename_profile[0],idx_lon,idx_lat))
    # idx=np.intersect1d(idx_lon, idx_lat)
    # print(idx)
    pressure_EcoPart_tmp=abs(pressure_EcoPart[sel])
    pressure_Coriolis_tmp=pres[i,:]
    ### Temperature interpolation
    temp_Coriolis_tmp=temp[i,:]
    selnan = temp_Coriolis_tmp != 99999
    fz=interp1d(pressure_Coriolis_tmp[selnan],temp_Coriolis_tmp[selnan],fill_value="extrapolate")#,kind='cubic')
    temp_EcoPart_tmp=fz(pressure_EcoPart_tmp)
    # fig = plt.figure(1, figsize=(12, 8))
    # plt.plot(temp_Coriolis_tmp,pressure_Coriolis_tmp , 'b')
    # plt.plot(temp_EcoPart_tmp,pressure_EcoPart_tmp,'r.')
    # plt.xlim(0,30);plt.ylim(0,500)
    temp_EcoPart[sel]=temp_EcoPart_tmp
    ### Dissolved oxygen interpolation
    doxy_Coriolis_tmp = doxy[i, :]
    selnan=doxy_Coriolis_tmp!=99999
    fz = interp1d(pressure_Coriolis_tmp[selnan], doxy_Coriolis_tmp[selnan],fill_value="extrapolate")  # ,kind='cubic')
    #sel_pressure_range=pressure_Coriolis_tmp[selnan].min()<=pressure_EcoPart_tmp<pressure_Coriolis_tmp[selnan].max()
    doxy_EcoPart_tmp = fz(pressure_EcoPart_tmp)
    doxy_EcoPart[sel] = doxy_EcoPart_tmp

data2append=np.ones((sel_filename.__len__(),2))*np.nan
data2append[sel_filename,:]=np.concatenate((temp_EcoPart.reshape(temp_EcoPart.size,1),doxy_EcoPart.reshape(doxy_EcoPart.size,1)),axis=1)
data2append=data2append.astype('float32')
column_labels=['Temperature [degrees Celsius]','Doxy [micromol/kg]']

#data=pd.read_csv(filename_ecopart, sep='\t', header=0)
for i in range(0,column_labels.__len__()):
    data[column_labels[i]]=data2append[:,i]

data.to_csv(filename_ecopart,sep='\t',index=False)

