import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import netCDF4 as nc
import seawater as sw
import gsw
import pickle
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from lin_fit import lin_fit
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
day0=datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime(2021,8,18)        # final date for the carbon budget calculation
ndays=(dayf-day0).days          # number of days
depth0=200                      # starting depth
depthf=600                      # final depth
layer_thickness=depthf-depth0   # thickness of the layer considered
delta_depth=15                  # around of the depth which I consider when extracting the flux
Oxy2C=0.89                      # to convert from mol of oxygen to mol of carbon
mol2gC=12.0107                  # to convert from mol of carbon to grams of carbon
day0_float=calendar.timegm(day0.timetuple())
dayf_float=calendar.timegm(dayf.timetuple())

########################################################################################################################
# I load and process data
########################################################################################################################

#I load the file with the flux and POC
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
Date_Time=np.array(data['Date_Time'][sel_filename])
depth=np.array(data['Depth [m]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])
bbp_POC=np.array(data['bbp POC [mgC/m3]'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I select the data only in the prescribed period
list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]

########################################################################################################################
# Here I calculate the integrated POC (i.e., MiP+MaP+bbp).
########################################################################################################################
MiP_POC_depth0_depthf=np.array([])
MaP_POC_depth0_depthf=np.array([])
bbp_POC_depth0_depthf=np.array([])
i=0
for i in range(0,list_dates.size):
    sel = Date_Num == list_dates[i]
    MiP_POC_tmp=MiP_POC[sel]
    MaP_POC_tmp=MaP_POC[sel]
    bbp_POC_tmp=MaP_POC[sel]
    y=depth[sel]
    sel2 = (~np.isnan(MiP_POC_tmp))&(~np.isnan(MaP_POC_tmp))&(~np.isnan(bbp_POC_tmp))
    MiP_POC_tmp = MiP_POC_tmp[sel2];MaP_POC_tmp = MaP_POC_tmp[sel2];bbp_POC_tmp = MaP_POC_tmp[sel2];y2 = y[sel2]
    #Here I select the values between depth0 and depthf
    sel_depth0_depthf = (np.abs(y) >= depth0) & (np.abs(y) < depthf)
    MiP_POC_depth0_depthf=np.append(MiP_POC_depth0_depthf,np.mean(MiP_POC_tmp[sel_depth0_depthf]))
    MaP_POC_depth0_depthf=np.append(MaP_POC_depth0_depthf,np.mean(MaP_POC_tmp[sel_depth0_depthf]))
    bbp_POC_depth0_depthf=np.append(bbp_POC_depth0_depthf,np.mean(bbp_POC_tmp[sel_depth0_depthf]))

Integrated_POC_mgC_m3 = MiP_POC_depth0_depthf + MaP_POC_depth0_depthf + bbp_POC_depth0_depthf
########################################################################################################################
# Here I extract the flux values at depth0 and depthf. To do so, (i) I filter it with a savgol function, then (ii) I
# interpolate it over a regular grid in time and space. This step is necessary to have the flux at 600 m, because some
# profiles only reach 400 m. Finally, (iii) I extract the flux values at depth0 and depthf
########################################################################################################################

##############################################
# Step 1 and 2, filter and interpolation
Flux_filtered=np.array([]);depth_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2) > 0:
        z = savgol_filter(z, 5, 1)
        Flux_filtered = np.concatenate((Flux_filtered, z))
        Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
        depth_Flux_filtered = np.concatenate((depth_Flux_filtered, y2))

# I define the x and y arrays for the Flux interpolation
x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
y_filtered = np.linspace(depth_Flux_filtered.min(), depth_Flux_filtered.max(), 99)
x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
# I interpolate
Flux_interp = griddata((Date_Num_Flux_filtered, depth_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")

##############################################
# Step 3, flux extraction at depth0 and depthf
Flux_depth0 = np.array([]);Flux_depthf = np.array([])

for i in range(0,Flux_interp.shape[1]):
    z=Flux_interp[:,i]
    sel_layer = (np.abs(y_filtered) >= depth0-delta_depth) & (np.abs(y_filtered) < depth0+delta_depth)
    Flux_depth0 = np.append(Flux_depth0, np.mean(z[sel_layer]))
    sel_layer = (np.abs(y_filtered) >= depthf-delta_depth) & (np.abs(y_filtered) < depthf+delta_depth)
    Flux_depthf = np.append(Flux_depthf, np.mean(z[sel_layer]))

########################################################################################################################
# Here I calculate the carbon consumption rate due to PARR
########################################################################################################################

############### I load Coriolis data with the oxygen information
filename='6903095_Sprof.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename))

lon=np.array(ds.variables['LONGITUDE'])
lat=np.array(ds.variables['LATITUDE'])
Date_Num=np.array(ds.variables['JULD'])
temp=np.array(ds.variables['TEMP_ADJUSTED'])
pres=np.array(ds.variables['PRES_ADJUSTED'])
psal=np.array(ds.variables['PSAL_ADJUSTED'])
doxy=np.array(ds.variables['DOXY_ADJUSTED'])

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

############### I load the PARR data from Ecopart
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

Date_Time_PARR=np.array(data['Date_Time'][sel_filename])
PARR_nmol_l_h=np.array(data['Respi_nmolO2_l_h'][sel_filename])
depth_PARR=np.array(data['Depth [m]'][sel_filename])
dens_PARR=np.array(data['Potential density [kg/m3]'][sel_filename])

# I convert the PARR measured in micromol/kg/day
PARR_micromol_kg_day=PARR_nmol_l_h.copy()/1000*24/(dens_PARR/1000)

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num_PARR=np.r_[0:Date_Time_PARR.size]
for i in range(0,Date_Time_PARR.size):
    date_time_obj = datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates_PARR=np.unique(Date_Num_PARR)

############### Here I start the loop on the different isopycnal values chosen for the study of the oxygen profile
delta_rho=0.05
reference_isopycnal_list=np.r_[1026.8:1027.50001:delta_rho]
Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+ndays]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))
depth_isopycnal=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
slopes_list_doxy=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
list_depth_PARR=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
list_dens_PARR=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))
PARR_list=np.squeeze(np.zeros((reference_isopycnal_list.size,1)))

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
    PARR_isopycnal=np.array([]);depth_PARR_tmp=np.array([]);dens_PARR_tmp=np.array([])
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

        # Here I proceed only if the date is inside the Date_Num_limit fixed
        if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

            # Here I extract the oxygen along the isopycnal
            sel_layer = (y >= reference_isopycnal_down) & (y < reference_isopycnal_up)
            if np.sum(sel_layer) > 0:  # If sel_layer has some True values, then I take as the doxy of this isopycnal the mean of the doxy values in correspondence with these layers
                doxy_tmp = np.mean(z[sel_layer])
                depth_isopycnal_tmp2 = np.mean(d[sel_layer])
                doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, depth_isopycnal_tmp2)
            else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the oxygen concentration at that depth
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
            if np.sum(sel_layer_PARR) > 0:  # If sel_layer_PARR has some True values, then I take as the PARR of this isopycnal the mean of the PARR values in correspondence with these layers
                PARR_tmp = np.mean(z_PARR[sel_layer_PARR])
                depth_PARR_tmp2 = np.mean(d_PARR[sel_layer_PARR])
                dens_PARR_tmp2 = np.mean(y_PARR[sel_layer_PARR])
                PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                depth_PARR_tmp = np.append(depth_PARR_tmp, depth_PARR_tmp2)
                dens_PARR_tmp = np.append(dens_PARR_tmp, dens_PARR_tmp2)
            else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the PARR at that depth
                depth_PARR_tmp2 = np.array([])
                dens_PARR_tmp2 = np.array([])
                for iy in range(0, y_PARR.size - 1):
                    if y_PARR[iy] <= reference_isopycnal < y_PARR[iy + 1]:
                        dist = (reference_isopycnal - y_PARR[iy]) / (y_PARR[iy + 1] - y_PARR[iy])
                        PARR_tmp = z_PARR[iy] + (z_PARR[iy + 1] - z_PARR[iy]) * dist
                        d_tmp = d_PARR[iy] + (d_PARR[iy + 1] - d_PARR[iy]) * dist
                        PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                        depth_PARR_tmp2 = np.append(depth_PARR_tmp2, d_tmp)
                        dens_PARR_tmp2 = np.append(dens_PARR_tmp2, reference_isopycnal)
                if depth_PARR_tmp2.size > 0:
                    depth_PARR_tmp = np.append(depth_PARR_tmp, np.mean(depth_PARR_tmp2))
                if dens_PARR_tmp2.size > 0:
                    dens_PARR_tmp = np.append(dens_PARR_tmp, np.mean(dens_PARR_tmp2))

    #Here I interpolate, so that the slope gives me the respiration rate (in micromol/kg per day)
    (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)
    depth_isopycnal[i_isop]=np.mean(depth_isopycnal_tmp)
    slopes_list_doxy[i_isop]=interpol.slope

    PARR_list[i_isop] = np.mean(PARR_isopycnal)
    list_depth_PARR[i_isop] = np.mean(depth_PARR_tmp)
    list_dens_PARR[i_isop] = np.mean(dens_PARR_tmp)

# I convert the PARR and the oxygen respiration rates (in micromolO2/kg/d) to the total amount of carbon consumption
# between depth0 and depthf, and between day0 and dayf (in mgC/m2)
# *Oxy2C -> to micromolC/kg/d
# *mol2gC -> to microgC/kg/d
# /1000 -> to mgC/kg/d
# *density -> to mgC/m3/d
# *layer_thickness*ndays -> to mgC/m2


POC_resp_mgC_m2=PARR_list*list_dens_PARR*Oxy2C*mol2gC*layer_thickness*ndays/1000
O2_resp_mgC_m2=slopes_list_doxy*reference_isopycnal_list*Oxy2C*mol2gC*layer_thickness*ndays/1000

########################################################################################################################
# Here I save the data
########################################################################################################################

dictionary_data = {"Integrated_POC_mgC_m3": Integrated_POC_mgC_m3, "list_dates_Integrated_POC": list_dates,
                   "Flux_depth0": Flux_depth0,"Flux_depthf": Flux_depthf, "Date_Num_Flux": x_filtered,
                   "POC_resp_mgC_m2": POC_resp_mgC_m2,"depth_POC_resp": list_depth_PARR, "dens_POC_resp": list_dens_PARR,
                   "O2_resp_mgC_m2": O2_resp_mgC_m2, "dens_02_resp": reference_isopycnal_list, "depth_02_resp": depth_isopycnal,
                   "day0": day0, "dayf": dayf,"ndays": ndays,"dayf_float": dayf_float,"day0_float": day0_float,
                   "depth0": depth0,"depthf": depthf,"layer_thickness": layer_thickness,
                   "delta_depth": delta_depth, "Oxy2C": Oxy2C,"mol2gC": mol2gC}

a_file = open("%s/an23/data_an23.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()

