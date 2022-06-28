import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4 as nc
import seawater as sw
import gsw
import sys
from scipy.interpolate import UnivariateSpline,griddata
import pickle
from pathlib import Path
home = str(Path.home())
#globals().clear()
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from lin_fit import lin_fit
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

# To update the data, please run:
# os.system("python Download_BGC_variables.py")
filename='6903095_Sprof_old.nc'

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

#Here I start the loop on the different isopycnal values chosen for the study of the oxygen profile
delta_rho=0.01
reference_isopycnal_list=np.r_[1026.8:1027.50001:delta_rho]
reference_isopycnal_list_plot=np.r_[0:reference_isopycnal_list.size+1:10]
Date_Num_limit=np.array([Date_Num.min(),Date_Num.max()]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))

depth_isopycnal=np.array([])
doxy_RR=np.array([]);Date_Num_RR=np.array([])
i_isop=0
for i_isop in range(0,reference_isopycnal_list.size):
    #Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
    #isopycnal, given by reference_isopycnal
    reference_isopycnal=reference_isopycnal_list[i_isop]
    reference_isopycnal_down=reference_isopycnal-delta_rho/2
    if reference_isopycnal_down<reference_isopycnal_list[0]:  reference_isopycnal_down=reference_isopycnal_list[0]
    reference_isopycnal_up=reference_isopycnal+delta_rho/2
    if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]
    parameter=doxy
    doxy_isopycnal=np.array([]);depth_isopycnal_tmp=np.array([]);Date_Num_isopycnal=np.array([])
    i=0
    for i in range(0,parameter.shape[0]):
        sel = (parameter[i, :] != 99999) & (dens[i, :] != 99999)
        z=parameter[i,sel];y=dens[i,sel];d=depth[i,sel]
        #fig = plt.figure(1);plt.plot(y,-d);plt.show();
        #plt.grid(color='k', linestyle='dashed', linewidth=0.5);plt.savefig('%s/GIT/AC_Agulhas_eddy_2021/Plots/an09/00_density_vs_depth_an09.pdf' % home);plt.close()
        if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:
            sel_layer = (y >= reference_isopycnal_down) & (y < reference_isopycnal_up)
            if np.sum(sel_layer)>0: #If sel_layer has some True values, then I take as the doxy of this isopycnal the mean of the doxy values in correspondence with these layers
                doxy_tmp = np.mean(z[sel_layer])
                depth_isopycnal_tmp2 = np.mean(d[sel_layer])
                doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, depth_isopycnal_tmp2)
            else: # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the oxygen concentration at that depth
                depth_isopycnal_tmp2 = np.array([])
                for iy in range(0, y.size - 1):
                    if y[iy] <= reference_isopycnal < y[iy + 1]:
                        dist = (reference_isopycnal - y[iy]) / (y[iy + 1] - y[iy])
                        doxy_tmp = z[iy] + (z[iy + 1] - z[iy]) * dist
                        d_tmp = d[iy] + (d[iy + 1] - d[iy]) * dist
                        doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                        depth_isopycnal_tmp2 = np.append(depth_isopycnal_tmp2, d_tmp)
                        Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                if depth_isopycnal_tmp2.size > 0:
                    depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, np.mean(depth_isopycnal_tmp2))

    if doxy_isopycnal.size>3:
        f=UnivariateSpline(Date_Num_isopycnal, doxy_isopycnal,k=3,s=round(Date_Num_isopycnal.size*1.5))
        # f1 = UnivariateSpline(Date_Num_isopycnal, doxy_isopycnal,k=3,s=round(Date_Num_isopycnal.size*1.2))
        # plt.scatter(Date_Num_isopycnal, doxy_isopycnal)
        # plt.plot(Date_Num_isopycnal, f(Date_Num_isopycnal), 'r')
        # plt.plot(Date_Num_isopycnal, f1(Date_Num_isopycnal), 'b')
        depth_isopycnal_tmp=np.squeeze(np.ones((1,Date_Num_isopycnal.size))*np.mean(depth_isopycnal_tmp))
        depth_isopycnal = np.append(depth_isopycnal, depth_isopycnal_tmp)
        fD =f.derivative()
        doxy_RR_tmp=fD(Date_Num_isopycnal)
        doxy_RR = np.append(doxy_RR, doxy_RR_tmp)
        Date_Num_RR = np.append(Date_Num_RR, Date_Num_isopycnal)

# I define the x and y arrays for the contourf plot
x_parameter = np.linspace(Date_Num_RR.min(), Date_Num_RR.max(), 100)
y1_parameter = np.linspace(depth_isopycnal.min(), depth_isopycnal.max(), 50)
# I interpolate
x_parameter_g, y_parameter_g = np.meshgrid(x_parameter, y1_parameter)
doxy_RR_interp_depth = griddata((Date_Num_RR, depth_isopycnal), doxy_RR,(x_parameter_g, y_parameter_g), method="nearest")

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = y1_parameter.min(),600
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num_RR.min(), Date_Num_RR.max()))
doxy_RR_plot=doxy_RR_interp_depth.copy()
doxy_RR_plot[doxy_RR_plot>0.0]=0.0
doxy_RR_plot[doxy_RR_plot<-0.07]=-0.07
ax_1 = plot2 = plt.contourf(x_parameter,y1_parameter, doxy_RR_plot,levels=20,cmap='Blues_r')#cmap='RdBu')#,vmin=-0.05)
plt.gca().invert_yaxis()
# draw colorbar
cbar = plt.colorbar(plot2)
cbar.ax.set_ylabel('Oxygen respiration rate ($\mu$mol kg$^{-1}$d$^{-1}$)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
#plt.title('%smm' % NP_sizeclass, fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num_RR.min(),Date_Num_RR.max(),nxticks)
xticklabels=[]
for i in xticks:
    date_time_obj = date_reference + datetime.timedelta(days=i)
    xticklabels.append(date_time_obj.strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an10/OxygenRespirationRate_vs_time_and_depth_an10.pdf',dpi=200)
plt.close()

# I save the data used for the plot
dictionary_data = {"x_parameter": x_parameter, "y1_parameter": y1_parameter, "doxy_RR_interp_depth": doxy_RR_interp_depth}
a_file = open("%s/an10/data_an10.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()