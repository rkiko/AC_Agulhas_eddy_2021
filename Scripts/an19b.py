import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4 as nc
import seawater as sw
import gsw
from pathlib import Path
from itertools import compress
import pandas as pd
import calendar
import sys
home = str(Path.home())
#globals().clear()
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from lin_fit import lin_fit
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts' % home) #changes directory
actualdir=os.getcwd()
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home



#######################################################################################################################
############### I load Coriolis data with the oxygen information
#######################################################################################################################

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
doxy=np.array(ds.variables['DOXY_ADJUSTED'])
doxy_qc=np.array(ds.variables['DOXY_ADJUSTED_QC'])
doxy_qc_profile=np.array(ds.variables['PROFILE_DOXY_QC'])

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
if np.sum(doxy==99999)==doxy.size:
    print('Taking non adjusted oxygen')
    doxy = np.array(ds.variables['DOXY'])
    doxy_qc = np.array(ds.variables['DOXY_QC'])

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



#######################################################################################################################
############### I load the PARR data from Ecopart
#######################################################################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_mip_map_flux_data.tsv' % home
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
lon_PARR=np.array(data['Longitude'][sel_filename])
lat_PARR=np.array(data['Latitude'][sel_filename])
Date_Time_PARR=np.array(data['Date_Time'][sel_filename])
PARR_nmol_l_h=np.array(data['Respi_nmolO2_l_h'][sel_filename])
depth_PARR=np.array(data['Depth [m]'][sel_filename])
dens_PARR=np.array(data['Potential density [kg/m3]'][sel_filename])
bbp_POC=np.array(data['bbp POC [mgC/m3]'][sel_filename])

#I convert the PARR measured in micromol/kg/day
PARR_micromol_kg_day=PARR_nmol_l_h.copy()/1000*24/(dens_PARR/1000)

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num_PARR=np.r_[0:Date_Time_PARR.size]
for i in range(0,Date_Time_PARR.size):
    date_time_obj = datetime.datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates_PARR=np.unique(Date_Num_PARR)
#######################################################################################################################
############### Here I start the loop on the different isopycnal values chosen for the study of the oxygen profile
#######################################################################################################################
if list_dates_PARR.size!=doxy.shape[0]: raise ValueError('The number of PARR profiles does not coincide with the number of oxygen profiles')

delta_rho=0.05
delta_rho_plot=0.1 #every how many delta_rho range I do a plot
reference_isopycnal_list=np.r_[1026.8:1027.50001:delta_rho]
reference_isopycnal_list_plot=np.r_[0:reference_isopycnal_list.size+1:np.round(delta_rho_plot/delta_rho)]
Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+71]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+71))

depth_isopycnal=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
slopes_list_doxy=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));slopes_ci_list_doxy=np.zeros((reference_isopycnal_list.size,2))
signif_list_doxy=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));signif_label_list_doxy=[]
R2_list_doxy=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));R_list_doxy=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))

list_depth_PARR=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
PARR_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
PARR_std_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))

slopes_list_bbp=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));slopes_ci_list_bbp=np.zeros((reference_isopycnal_list.size,2))
signif_list_bbp=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));signif_label_list_bbp=[]
R2_list_bbp=np.squeeze(np.zeros((reference_isopycnal_list.size,1)));R_list_bbp=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))

i_isop=0
for i_isop in range(0,reference_isopycnal_list.size):
    #Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
    #isopycnal, given by reference_isopycnal, and then the PARR concentration at that isopycnal
    reference_isopycnal=reference_isopycnal_list[i_isop]
    reference_isopycnal_down=reference_isopycnal-delta_rho/2
    if reference_isopycnal_down<reference_isopycnal_list[0]:  reference_isopycnal_down=reference_isopycnal_list[0]
    reference_isopycnal_up=reference_isopycnal+delta_rho/2
    if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]
    doxy_isopycnal=np.array([]);depth_isopycnal_tmp=np.array([]);Date_Num_isopycnal=np.array([])
    PARR_isopycnal=np.array([]);depth_PARR_tmp=np.array([])
    bbp_isopycnal=np.array([]);depth_bbp_tmp=np.array([]);Date_Num_isopycnal_bbp=np.array([])
    i=0
    for i in range(0,doxy.shape[0]):
        #Here, for the i-th profile, I select the oxygen, density and depth profiles of Coriolis data, excluding the nan values
        sel = (doxy[i, :] != 99999) & (dens[i, :] != 99999)
        z=doxy[i,sel];y=dens[i,sel];d=depth[i,sel]
        # Here, for the i-th profile, I select the PARR, density and depth profiles of Ecopart data, excluding the nan values. I do the same for the bbp
        sel_PARR=Date_Num_PARR==list_dates_PARR[i]
        z_PARR=PARR_micromol_kg_day[sel_PARR];y_PARR=dens_PARR[sel_PARR];d_PARR=depth_PARR[sel_PARR]
        z_bbp=bbp_POC[sel_PARR];y_bbp=dens_PARR[sel_PARR];d_bbp=depth_PARR[sel_PARR]
        sel_PARR = (~np.isnan(z_PARR)) & (~np.isnan(y_PARR))
        z_PARR = z_PARR[sel_PARR];y_PARR = y_PARR[sel_PARR];d_PARR = d_PARR[sel_PARR]
        sel_bbp = (~np.isnan(z_bbp)) & (~np.isnan(y_bbp))
        z_bbp = z_bbp[sel_bbp];y_bbp = y_bbp[sel_bbp];d_bbp = d_bbp[sel_bbp]
        # Here I proceed only if the date is inside the Date_Num_limit fixed
        if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

            # Here I extract the oxygen along the isopycnal
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

            # Here I extract the PARR along the isopycnal
            sel_layer_PARR = (y_PARR >= reference_isopycnal_down) & (y_PARR < reference_isopycnal_up)
            if np.sum(sel_layer_PARR)>0: #If sel_layer_PARR has some True values, then I take as the PARR of this isopycnal the mean of the PARR values in correspondence with these layers
                PARR_tmp = np.mean(z_PARR[sel_layer_PARR])
                depth_PARR_tmp2 = np.mean(d_PARR[sel_layer_PARR])
                PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                depth_PARR_tmp = np.append(depth_PARR_tmp, depth_PARR_tmp2)
            else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the PARR at that depth
                depth_PARR_tmp2 = np.array([])
                for iy in range(0, y_PARR.size - 1):
                    if y_PARR[iy] <= reference_isopycnal < y_PARR[iy + 1]:
                        dist = (reference_isopycnal - y_PARR[iy]) / (y_PARR[iy + 1] - y_PARR[iy])
                        PARR_tmp = z_PARR[iy] + (z_PARR[iy + 1] - z_PARR[iy]) * dist
                        d_tmp = d_PARR[iy] + (d_PARR[iy + 1] - d_PARR[iy]) * dist
                        PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                        depth_PARR_tmp2 = np.append(depth_PARR_tmp2, d_tmp)
                if depth_PARR_tmp2.size > 0:
                    depth_PARR_tmp = np.append(depth_PARR_tmp, np.mean(depth_PARR_tmp2))

            # Here I extract the bbp_POC along the isopycnal
            sel_layer_bbp = (y_bbp >= reference_isopycnal_down) & (y_bbp < reference_isopycnal_up)
            if np.sum(sel_layer_bbp)>0: #If sel_layer_bbp has some True values, then I take as the bbp of this isopycnal the mean of the bbp values in correspondence with these layers
                bbp_tmp = np.mean(z_bbp[sel_layer_bbp])
                depth_bbp_tmp2 = np.mean(d_bbp[sel_layer_bbp])
                bbp_isopycnal = np.append(bbp_isopycnal, bbp_tmp)
                Date_Num_isopycnal_bbp = np.append(Date_Num_isopycnal_bbp, Date_Num[i])
                depth_bbp_tmp = np.append(depth_bbp_tmp, depth_bbp_tmp2)
            else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the bbp at that depth
                depth_bbp_tmp2 = np.array([])
                for iy in range(0, y_bbp.size - 1):
                    if y_bbp[iy] <= reference_isopycnal < y_bbp[iy + 1]:
                        dist = (reference_isopycnal - y_bbp[iy]) / (y_bbp[iy + 1] - y_bbp[iy])
                        bbp_tmp = z_bbp[iy] + (z_bbp[iy + 1] - z_bbp[iy]) * dist
                        d_tmp = d_bbp[iy] + (d_bbp[iy + 1] - d_bbp[iy]) * dist
                        bbp_isopycnal = np.append(bbp_isopycnal, bbp_tmp)
                        depth_bbp_tmp2 = np.append(depth_bbp_tmp2, d_tmp)
                        Date_Num_isopycnal_bbp = np.append(Date_Num_isopycnal_bbp, Date_Num[i])
                if depth_bbp_tmp2.size > 0:
                    depth_bbp_tmp = np.append(depth_bbp_tmp, np.mean(depth_bbp_tmp2))

    #Here I interpolate, so that the slope gives me the respiration rate (in micromol/kg per day)
    (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)

    depth_isopycnal[i_isop]=np.mean(depth_isopycnal_tmp)
    if depth_isopycnal[0]==3: print(np.mean(depth_isopycnal_tmp),i_isop)
    slopes_list_doxy[i_isop]=interpol.slope
    slopes_ci_list_doxy[i_isop,:]=slpe_ci
    signif_list_doxy[i_isop]=signif
    signif_label_list_doxy.append(signif_label)
    R2_list_doxy[i_isop] = interpol.rvalue ** 2
    R_list_doxy[i_isop] = interpol.rvalue

    # Here I interpolate, so that the slope gives me the carbon (i.e. bbp_POC) consumption rate  (in mmolC/kg per day)
    (interpol, slpe_ci, _, signif, signif_label) = lin_fit(Date_Num_isopycnal_bbp, bbp_isopycnal)

    slopes_list_bbp[i_isop] = interpol.slope
    slopes_ci_list_bbp[i_isop, :] = slpe_ci
    signif_list_bbp[i_isop] = signif
    signif_label_list_bbp.append(signif_label)
    R2_list_bbp[i_isop] = interpol.rvalue ** 2
    R_list_bbp[i_isop] = interpol.rvalue

    PARR_list[i_isop] = np.mean(PARR_isopycnal)
    PARR_std_list[i_isop] = np.std(PARR_isopycnal)
    list_depth_PARR[i_isop] = np.mean(depth_PARR_tmp)

    ##############################################################################
    ##############################################################################
    ######Plot

    # if np.sum(np.isin(reference_isopycnal_list_plot,i_isop))==0:    continue
    #Here I plot
    width, height = 0.8, 0.7
    fig = plt.figure(1, figsize=(12, 8))
    ax = fig.add_axes([0.12, 0.2, width, height])
    plot3 = plt.scatter(Date_Num_isopycnal_bbp, bbp_isopycnal, c='black')
    plot4 = plt.plot(np.linspace(Date_Num_limit[0], Date_Num_limit[1], 20),
                     np.linspace(Date_Num_limit[0] * interpol.slope + interpol.intercept,
                                 Date_Num_limit[1] * interpol.slope + interpol.intercept, 20), c='black')
    plt.text(Date_Num_limit[0]+(Date_Num_limit[1]-Date_Num_limit[0])*0.1,bbp_isopycnal.min()+(bbp_isopycnal.max()-bbp_isopycnal.min())*0.2,'R$^2$=%.2f, p=%0.2e' % (interpol.rvalue**2,interpol.pvalue), fontsize=16)
    plt.ylabel('bbp POC (mg/m$^3$)', fontsize=18)
    plt.title('bbp time series along isop. %.1f kg/m$^3$ (%d m), slope: %0.3f mg/m$^3$/d fit signif: %s' % (reference_isopycnal,depth_isopycnal[i_isop],interpol.slope / 0.89,signif_label), fontsize=15)
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
    plt.savefig('../Plots/an19b/02_bbpTimeSeries_isopycnal%0.2f_deltarho_%0.2f_an19b.pdf' % (reference_isopycnal,delta_rho), dpi=200)
    plt.close()
    ##############################################################################
    ##############################################################################

# I divide by 0.89 to convert the carbon (i.e. bbp_POC) consumption rate to oxygen consumption rate
slopes_list_bbp = slopes_list_bbp / 0.89
slopes_ci_list_bbp = slopes_ci_list_bbp / 0.89


#Here I plot the respiration rate vs the depth
fig, ax = plt.subplots(1,2,figsize=(6,8), gridspec_kw={'width_ratios': [4, 1]},sharey=True)
plt.subplots_adjust(wspace=0.2)
#ax = fig.add_axes([0.12, 0.2, width, height])
ax[0].plot(slopes_list_doxy, depth_isopycnal, 'k')#, linestyle='dashed')
ax[0].scatter(slopes_list_doxy, depth_isopycnal, c='black',s=5)
ax[0].fill_betweenx(depth_isopycnal,slopes_ci_list_doxy[:,1],slopes_ci_list_doxy[:,0],facecolor='b',color='gray',alpha=0.5,label='O$_2$')
ax[0].plot(-PARR_list, list_depth_PARR, 'b')#, linestyle='dashed')
ax[0].scatter(-PARR_list, list_depth_PARR, c='b',s=5)
ax[0].fill_betweenx(list_depth_PARR,-(PARR_list-PARR_std_list),-(PARR_list+PARR_std_list),facecolor='b',color='b',alpha=0.5,label='PARR')
ax[0].plot(slopes_list_bbp, depth_isopycnal, 'green')#, linestyle='dashed')
ax[0].scatter(slopes_list_bbp, depth_isopycnal, c='black',s=5)
ax[0].fill_betweenx(depth_isopycnal,slopes_ci_list_bbp[:,1],slopes_ci_list_bbp[:,0],facecolor='b',color='green',alpha=0.5,label='bbp O$_2$')
ax[0].legend(fontsize=12,loc='center left')
ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
ax[0].set_ylim(depth_isopycnal.min(),600)
ax[0].invert_yaxis()
ax[0].set_xticks(np.r_[-0.1:0.101:0.02])
ax[0].set_yticks(np.r_[200:601:50])
ax[0].set_xlabel('Respiration rate ($\mu$mol kg$^{-1}$d$^{-1}$)', fontsize=12)
ax[0].set_ylabel('Depth (m)', fontsize=12)
ax[0].tick_params(labelsize=12)
#set_xticklabels(fontsize=18),plt.yticks(fontsize=18)
ax[0].set_xlim(-0.08,0)
ax[0].set_title('Resp. rate vs depth (2021-04-13 to 2021-06-23)', fontsize=10)
ax[1].scatter(R2_list_doxy[R_list_doxy<0], depth_isopycnal[R_list_doxy<0], c='r',s=5)#, linestyle='dashed')
ax[1].set_xlim(0,1.3)
ax[1].axvline(1, linestyle='-', color='k')
ax[1].set_xticks(np.array([0,0.5,1.]))
ax[1].set_xticks(np.array([0.25,0.75]),minor=True)
ax[1].set_yticks(np.r_[150:601:50])
ax[1].grid(color='k', linestyle='dashed', linewidth=0.5,which='both')
ax[1].set_ylim(depth_isopycnal.min(),600)
ax[1].invert_yaxis()
ax[1].set_xlabel('R$^2$', fontsize=12)
ax[1].set_title('O$_2$ Fit significance', fontsize=10)
sel=(R_list_doxy<0) & (depth_isopycnal< 600)
xtext=np.squeeze(np.ones((1,depth_isopycnal[sel].size)))
ytext=depth_isopycnal[sel]
text_s=list(compress(signif_label_list_doxy, sel))
for i,j,k in zip(xtext,ytext,text_s):
    ax[1].text(i,j,k,va='top')
fig.savefig('../Plots/an19b/01_OxygenRespirationRate_vs_depth_deltarho_%0.2f_an19b.pdf' % delta_rho, dpi=200)
plt.close()

