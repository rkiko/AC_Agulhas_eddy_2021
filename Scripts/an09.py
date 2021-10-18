import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4 as nc
import seawater as sw
import gsw
from pathlib import Path
from itertools import compress
import sys
home = str(Path.home())
#globals().clear()
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from lin_fit import lin_fit
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

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
delta_rho_plot=0.1 #every how many delta_rho range I do a plot
reference_isopycnal_list=np.r_[1026.8:1027.50001:delta_rho]
reference_isopycnal_list_plot=np.r_[0:reference_isopycnal_list.size+1:np.round(delta_rho_plot/delta_rho)]
Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+127]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))

depth_isopycnal=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
slopes_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));slopes_ci_list=np.zeros((reference_isopycnal_list.size,2))
signif_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));signif_label_list=[]
R2_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));R_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
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
    for i in range(0,doxy.shape[0]):
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

    #Here I interpolate, so that the slope gives me the respiration rate (in micromol/kg per day)
    (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)

    depth_isopycnal[i_isop]=np.mean(depth_isopycnal_tmp)
    if depth_isopycnal[0]==3: print(np.mean(depth_isopycnal_tmp),i_isop)
    slopes_list[i_isop]=interpol.slope
    slopes_ci_list[i_isop,:]=slpe_ci
    signif_list[i_isop]=signif
    signif_label_list.append(signif_label)
    R2_list[i_isop] = interpol.rvalue ** 2
    R_list[i_isop] = interpol.rvalue

    if np.sum(np.isin(reference_isopycnal_list_plot,i_isop))==0:    continue
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
    plt.savefig('../Plots/an09/02_OxygenTimeSeries_isopycnal%0.1f_deltarho_%0.2f_an09.pdf' % (reference_isopycnal,delta_rho), dpi=200)
    plt.close()

#Here I plot the respiration rate vs the depth
width, height = 0.8, 0.7
fig, ax = plt.subplots(1,2,figsize=(12,8), gridspec_kw={'width_ratios': [8, 1]},sharey=True)
plt.subplots_adjust(wspace=0.1)
#ax = fig.add_axes([0.12, 0.2, width, height])
ax[0].plot(slopes_list, depth_isopycnal, 'r')#, linestyle='dashed')
ax[0].scatter(slopes_list, depth_isopycnal, c='black',s=5)
ax[0].axvline(0, linestyle='-', color='k')
ax[0].fill_betweenx(depth_isopycnal,slopes_ci_list[:,0],slopes_ci_list[:,1],facecolor='b',color='gray',alpha=0.5)
ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
ax[0].set_ylim(depth_isopycnal.min(),600)
ax[0].invert_yaxis()
ax[0].set_xticks(np.r_[-0.1:0.101:0.02])
ax[0].set_yticks(np.r_[150:601:50])
ax[0].set_xlabel('Oxygen respiration rate ($\mu$mol kg$^{-1}$d$^{-1}$)', fontsize=18)
ax[0].set_ylabel('Depth (m)', fontsize=18)
ax[0].tick_params(labelsize=18)
#set_xticklabels(fontsize=18),plt.yticks(fontsize=18)
ax[0].set_xlim(-0.1,0.1)
ax[0].set_title('O$_2$ resp. rate vs depth (2021-04-13 to 2021-08-18)', fontsize=18)
ax[1].scatter(R2_list[R_list<0], depth_isopycnal[R_list<0], c='r',s=3)#, linestyle='dashed')
ax[1].set_xlim(0,1.3)
ax[1].axvline(1, linestyle='-', color='k')
ax[1].set_xticks(np.array([0,0.5,1.]))
ax[1].set_xticks(np.array([0.25,0.75]),minor=True)
ax[1].set_yticks(np.r_[150:601:50])
ax[1].grid(color='k', linestyle='dashed', linewidth=0.5,which='both')
ax[1].set_ylim(depth_isopycnal.min(),600)
ax[1].invert_yaxis()
ax[1].set_xlabel('R$^2$', fontsize=18)
ax[1].set_title('Fit significance', fontsize=12)
sel=(R_list<0) & (depth_isopycnal< 600)
xtext=np.squeeze(np.ones((1,depth_isopycnal[sel].size)))
ytext=depth_isopycnal[sel]
text_s=list(compress(signif_label_list, sel))
for i,j,k in zip(xtext,ytext,text_s):
    ax[1].text(i,j,k,va='top')
fig.savefig('../Plots/an09/01_OxygenRespirationRate_vs_depth_deltarho_%0.2f_an09.pdf' % delta_rho, dpi=200)
plt.close()

