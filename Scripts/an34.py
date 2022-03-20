import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime,timedelta
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from scipy.stats import norm
import netCDF4 as nc
import seawater as sw
import gsw
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
from lin_fit import lin_fit
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
day0=datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime(2021,9,24)        # starting date for the carbon budget calculation
dayIM1=datetime(2021,6,20)      # Start intermediate date
dayIM2=datetime(2021,8,18)      # Final intermediate date
ndays1=(dayIM1-day0).days       # number of days
ndays2=(dayIM2-day0).days       # number of days
delta_days=5                    # every how many days I put the intermediate date for the carbon budget calculation
depth00=200                     # starting depth
layer_thickness=100             # thickness of the layer considered
delta_depth=25                  # every time I do a loop, how much I do increase depth0
depthff=600                     # final maximal depth investigated

depth0_list=np.r_[depth00:depthff-layer_thickness+0.1:delta_depth]

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Here I define the function
########################################################################################################################
########################################################################################################################
########################################################################################################################

def carbon_budget_calculation(depth0,depthf,day0,dayf):
    ########################################################################################################################
    # Starting parameters
    ########################################################################################################################
    # dayf = day0+timedelta(days=ndays) # final date for the carbon budget calculation
    ndays = (dayf - day0).days  # number of days
    delta_depth_flux = 15  # around of the depth which I consider when extracting the flux
    Oxy2C = 0.89  # to convert from mol of oxygen to mol of carbon
    mol2gC = 12.0107  # to convert from mol of carbon to grams of carbon
    day0_float = calendar.timegm(day0.timetuple())
    dayf_float = calendar.timegm(dayf.timetuple())

    ########################################################################################################################
    # I load and process data
    ########################################################################################################################

    #I load the file with the flux and POC
    filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_356.tsv' % home
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
    Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
    Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
    Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
    MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
    MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
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
    n_profiles=list_dates.size

    ########################################################################################################################
    # Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid in time and space. This step is necessary to have MiP+MaP+bbp at 600 m, because some
    # profiles only reach 400 m. Finally, (iii) I extract the mean MiP+MaP+bbp values between depth0 and depthf and between
    # day0 and dayf (I obtain a time series)
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    MiP_filtered=np.array([]);depth_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
    MiP_extended_filtered=np.array([]);depth_MiP_extended_filtered=np.array([]);Date_Num_MiP_extended_filtered=np.array([])
    MaP_filtered=np.array([]);depth_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
    bbp_filtered=np.array([]);depth_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])

    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=MiP_POC[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_filtered = np.concatenate((MiP_filtered, z))
            Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x2))
            depth_MiP_filtered = np.concatenate((depth_MiP_filtered, y2))
        z=MiP_POC_extended[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_extended_filtered = np.concatenate((MiP_extended_filtered, z))
            Date_Num_MiP_extended_filtered = np.concatenate((Date_Num_MiP_extended_filtered, x2))
            depth_MiP_extended_filtered = np.concatenate((depth_MiP_extended_filtered, y2))
        z=MaP_POC[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MaP_filtered = np.concatenate((MaP_filtered, z))
            Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x2))
            depth_MaP_filtered = np.concatenate((depth_MaP_filtered, y2))
        z=bbp_POC[sel];x=Date_Num[sel];y=depth[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            bbp_filtered = np.concatenate((bbp_filtered, z))
            Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, x2))
            depth_bbp_filtered = np.concatenate((depth_bbp_filtered, y2))

    # I define the x and y arrays for the MiP+MaP+bbp interpolation
    x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), ndays)
    y_filtered = np.linspace(depth_bbp_filtered.min(), depth_bbp_filtered.max(), 100)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    MiP_interp = griddata((Date_Num_MiP_filtered, depth_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MiP_extended_interp = griddata((Date_Num_MiP_extended_filtered, depth_MiP_extended_filtered), MiP_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MaP_interp = griddata((Date_Num_MaP_filtered, depth_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_interp = griddata((Date_Num_bbp_filtered, depth_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")

    ##############################################
    # Step 3, I calculate the mean MiP+MaP+bbp (and std) between depth0 and depthf between day0 and dayf
    sel_depth0_depthf = (np.abs(y_filtered) >= depth0) & (np.abs(y_filtered) < depthf)
    MiP_POC_depth0_depthf=np.mean(MiP_interp[sel_depth0_depthf,:],0)
    MiP_POC_extended_depth0_depthf=np.mean(MiP_extended_interp[sel_depth0_depthf,:],0)
    MaP_POC_depth0_depthf=np.mean(MaP_interp[sel_depth0_depthf,:],0)
    bbp_POC_depth0_depthf=np.mean(bbp_interp[sel_depth0_depthf,:],0)

    MiP_POC_depth0_depthf_std = np.std(MiP_interp[sel_depth0_depthf, :], 0)
    MiP_POC_extended_depth0_depthf_std = np.std(MiP_extended_interp[sel_depth0_depthf, :], 0)
    MaP_POC_depth0_depthf_std = np.std(MaP_interp[sel_depth0_depthf, :], 0)
    bbp_POC_depth0_depthf_std = np.std(bbp_interp[sel_depth0_depthf, :], 0)

    Integrated_POC_mgC_m3 = MiP_POC_depth0_depthf + MaP_POC_depth0_depthf + bbp_POC_depth0_depthf
    Integrated_POC_extended_mgC_m3 = MiP_POC_extended_depth0_depthf + MaP_POC_depth0_depthf + bbp_POC_depth0_depthf
    Integrated_POC_mgC_m3_std = np.sqrt( MiP_POC_depth0_depthf_std**2 + MaP_POC_depth0_depthf_std**2 + bbp_POC_depth0_depthf_std**2 )
    Integrated_POC_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_depth0_depthf_std**2 + MaP_POC_depth0_depthf_std**2 + bbp_POC_depth0_depthf_std**2 )
    list_dates_Integrated_POC = x_filtered.copy()
    ########################################################################################################################
    # Here I extract the flux values at depth0 and depthf. To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid in time and space. This step is necessary to have the flux at 600 m, because some
    # profiles only reach 400 m. Finally, (iii) I extract the flux values at depth0 and depthf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    Flux_filtered=np.array([]);depth_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
    Flux_eta_b_filtered=np.array([]);depth_Flux_eta_b_filtered=np.array([]);Date_Num_Flux_eta_b_filtered=np.array([])
    Flux_extended_filtered=np.array([]);depth_Flux_extended_filtered=np.array([]);Date_Num_Flux_extended_filtered=np.array([])
    Flux_extended_filtered_eta_b=np.array([]);depth_Flux_extended_filtered_eta_b=np.array([]);Date_Num_Flux_extended_filtered_eta_b=np.array([])
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
        z=Flux_eta_b[sel];x=Date_Num[sel];y=depth[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_eta_b_filtered = np.concatenate((Flux_eta_b_filtered, z))
            Date_Num_Flux_eta_b_filtered = np.concatenate((Date_Num_Flux_eta_b_filtered, x2))
            depth_Flux_eta_b_filtered = np.concatenate((depth_Flux_eta_b_filtered, y2))
        z=Flux_extended[sel];x=Date_Num[sel];y=depth[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered = np.concatenate((Flux_extended_filtered, z))
            Date_Num_Flux_extended_filtered = np.concatenate((Date_Num_Flux_extended_filtered, x2))
            depth_Flux_extended_filtered = np.concatenate((depth_Flux_extended_filtered, y2))
        z=Flux_extended_eta_b[sel];x=Date_Num[sel];y=depth[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered_eta_b = np.concatenate((Flux_extended_filtered_eta_b, z))
            Date_Num_Flux_extended_filtered_eta_b = np.concatenate((Date_Num_Flux_extended_filtered_eta_b, x2))
            depth_Flux_extended_filtered_eta_b = np.concatenate((depth_Flux_extended_filtered_eta_b, y2))

    # I define the x and y arrays for the Flux interpolation
    x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
    y_filtered = np.linspace(depth_Flux_filtered.min(), depth_Flux_filtered.max(), 99)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    Flux_interp = griddata((Date_Num_Flux_filtered, depth_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_eta_b_interp = griddata((Date_Num_Flux_eta_b_filtered, depth_Flux_eta_b_filtered), Flux_eta_b_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_interp = griddata((Date_Num_Flux_extended_filtered, depth_Flux_extended_filtered), Flux_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_eta_b_interp = griddata((Date_Num_Flux_extended_filtered_eta_b, depth_Flux_extended_filtered_eta_b), Flux_extended_filtered_eta_b,(x_filtered_g, y_filtered_g), method="nearest")

    ##############################################
    # Step 3, flux extraction at depth0 and depthf
    Flux_depth0 = np.array([]);Flux_depthf = np.array([])
    Flux_eta_b_depth0 = np.array([]);Flux_eta_b_depthf = np.array([])
    Flux_extended_depth0 = np.array([]);Flux_extended_depthf = np.array([])
    Flux_extended_eta_b_depth0 = np.array([]);Flux_extended_eta_b_depthf = np.array([])

    for i in range(0,Flux_interp.shape[1]):
        z=Flux_interp[:,i]
        sel_layer = (np.abs(y_filtered) >= depth0-delta_depth_flux) & (np.abs(y_filtered) < depth0+delta_depth_flux)
        Flux_depth0 = np.append(Flux_depth0, np.mean(z[sel_layer]))
        sel_layer = (np.abs(y_filtered) >= depthf - delta_depth_flux) & (np.abs(y_filtered) < depthf + delta_depth_flux)
        Flux_depthf = np.append(Flux_depthf, np.mean(z[sel_layer]))

        z=Flux_eta_b_interp[:,i]
        sel_layer = (np.abs(y_filtered) >= depth0-delta_depth_flux) & (np.abs(y_filtered) < depth0+delta_depth_flux)
        Flux_eta_b_depth0 = np.append(Flux_eta_b_depth0, np.mean(z[sel_layer]))
        sel_layer = (np.abs(y_filtered) >= depthf - delta_depth_flux) & (np.abs(y_filtered) < depthf + delta_depth_flux)
        Flux_eta_b_depthf = np.append(Flux_eta_b_depthf, np.mean(z[sel_layer]))

        z=Flux_extended_interp[:,i]
        sel_layer = (np.abs(y_filtered) >= depth0 - delta_depth_flux) & (np.abs(y_filtered) < depth0 + delta_depth_flux)
        Flux_extended_depth0 = np.append(Flux_extended_depth0, np.mean(z[sel_layer]))
        sel_layer = (np.abs(y_filtered) >= depthf - delta_depth_flux) & (np.abs(y_filtered) < depthf + delta_depth_flux)
        Flux_extended_depthf = np.append(Flux_extended_depthf, np.mean(z[sel_layer]))

        z=Flux_extended_eta_b_interp[:,i]
        sel_layer = (np.abs(y_filtered) >= depth0 - delta_depth_flux) & (np.abs(y_filtered) < depth0 + delta_depth_flux)
        Flux_extended_eta_b_depth0 = np.append(Flux_extended_eta_b_depth0, np.mean(z[sel_layer]))
        sel_layer = (np.abs(y_filtered) >= depthf - delta_depth_flux) & (np.abs(y_filtered) < depthf + delta_depth_flux)
        Flux_extended_eta_b_depthf = np.append(Flux_extended_eta_b_depthf, np.mean(z[sel_layer]))


    ########################################################################################################################
    # Here I calculate the carbon consumption rate due to (i) oxygen consumption and (ii) PARR
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

    ############### I load the PARR correlated data (density, depth, and date) from Ecopart
    Date_Time_PARR=np.array(data['Date_Time'][sel_filename])
    depth_PARR=np.array(data['Depth [m]'][sel_filename])
    dens_PARR=np.array(data['Potential density [kg/m3]'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_PARR=np.r_[0:Date_Time_PARR.size]
    for i in range(0,Date_Time_PARR.size):
        date_time_obj = datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    list_dates_PARR=np.unique(Date_Num_PARR)

    #############################################
    ############### Loop on the different isopycnal values chosen for the study of the oxygen profile
    delta_rho=0.05
    reference_isopycnal_list=np.r_[1026.8:1027.50001:delta_rho]
    Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+ndays]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))
    depth_isopycnal = np.array([])
    slopes_list_doxy = np.array([])
    slopes_ci_list_doxy = np.array([])
    reference_isopycnal_list_sel = np.array([])

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
        i=0
        for i in range(0,doxy.shape[0]):
            #Here, for the i-th profile, I select the oxygen, density and depth profiles of Coriolis data, excluding the nan values
            sel = (doxy[i, :] != 99999) & (dens[i, :] != 99999)
            z=doxy[i,sel];y=dens[i,sel];d=depth[i,sel]

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

        # If I have at least three points I interpolate them, so that the slope gives me the respiration rate
        # (in micromol/kg per day)
        if Date_Num_isopycnal.size>2:
            (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)
            # fig = plt.figure(1, figsize=(12, 8))
            # plot1 = plt.scatter(Date_Num_isopycnal,doxy_isopycnal)
            depth_isopycnal = np.append(depth_isopycnal ,np.mean(depth_isopycnal_tmp))
            slopes_list_doxy = np.append(slopes_list_doxy ,interpol.slope)
            reference_isopycnal_list_sel = np.append(reference_isopycnal_list_sel ,reference_isopycnal)
            if i_isop==0:
                slopes_ci_list_doxy = np.reshape(slpe_ci.copy(),(1,2))
            else:
                slopes_ci_list_doxy = np.concatenate( ( slopes_ci_list_doxy ,np.reshape(slpe_ci,(1,2)) ), axis=0)

    O2_resp_mgC_m2_vs_depth=slopes_list_doxy.copy()*reference_isopycnal_list_sel*Oxy2C*mol2gC*layer_thickness*ndays/1000
    O2_resp_mgC_m2_vs_depth_ci=slopes_ci_list_doxy.copy()*( np.tile(reference_isopycnal_list_sel,(2,1)).T) *Oxy2C*mol2gC*layer_thickness*ndays/1000

    #############################################
    ############### Loop on different respiration types used to estimate PARR
    #List of the different Respiration types present in data
    list_Respi_types = [match for match in data.columns if "Respi" in match]
    nRespi= len(list_Respi_types)  #number of respiration types

    POC_resp_mgC_m2_vs_depth_list=np.zeros((reference_isopycnal_list.size,nRespi))
    POC_resp_mgC_m2_std_vs_depth_list = np.zeros((reference_isopycnal_list.size, nRespi))

    iRespi=0
    for iRespi in range(0,nRespi):
        PARR_nmol_l_h = np.array(data[list_Respi_types[iRespi]][sel_filename])
        # I convert the PARR measured in micromol/kg/day
        PARR_micromol_kg_day = PARR_nmol_l_h.copy() / 1000 * 24 / (dens_PARR / 1000)

        #############################################
        ############### Loop on the different isopycnal values chosen for the study of the PARR profile
        list_depth_PARR = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        list_dens_PARR = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        PARR_list = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))
        PARR_std_list = np.squeeze(np.zeros((reference_isopycnal_list.size, 1)))

        i_isop=0
        for i_isop in range(0,reference_isopycnal_list.size):
            #Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration at a given
            #isopycnal, given by reference_isopycnal, and then the PARR concentration at that isopycnal
            reference_isopycnal=reference_isopycnal_list[i_isop]
            reference_isopycnal_down=reference_isopycnal-delta_rho/2
            if reference_isopycnal_down<reference_isopycnal_list[0]:  reference_isopycnal_down=reference_isopycnal_list[0]
            reference_isopycnal_up=reference_isopycnal+delta_rho/2
            if reference_isopycnal_up > reference_isopycnal_list[-1]:  reference_isopycnal_up = reference_isopycnal_list[-1]
            PARR_isopycnal=np.array([]);depth_PARR_tmp=np.array([]);dens_PARR_tmp=np.array([])
            i=0
            for i in range(0,doxy.shape[0]):
                # Here, for the i-th profile, I select the PARR, density and depth profiles of Ecopart data, excluding the nan values. I do the same for the bbp
                sel_PARR=Date_Num_PARR==list_dates_PARR[i]
                z_PARR=PARR_micromol_kg_day[sel_PARR];y_PARR=dens_PARR[sel_PARR];d_PARR=depth_PARR[sel_PARR]
                sel_PARR = (~np.isnan(z_PARR)) & (~np.isnan(y_PARR))
                z_PARR = z_PARR[sel_PARR];y_PARR = y_PARR[sel_PARR];d_PARR = d_PARR[sel_PARR]

                # Here I proceed only if the date is inside the Date_Num_limit fixed
                if Date_Num_limit[0] <= Date_Num[i] <= Date_Num_limit[1]:

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

            PARR_list[i_isop] = np.mean(PARR_isopycnal)
            PARR_std_list[i_isop] = np.std(PARR_isopycnal)
            list_depth_PARR[i_isop] = np.mean(depth_PARR_tmp)
            list_dens_PARR[i_isop] = np.mean(dens_PARR_tmp)

        # I convert the PARR and the oxygen respiration rates (in micromolO2/kg/d) to the total amount of carbon consumption
        # between depth0 and depthf, and between day0 and dayf (in mgC/m2)
        # *Oxy2C -> to micromolC/kg/d
        # *mol2gC -> to microgC/kg/d
        # /1000 -> to mgC/kg/d
        # *density -> to mgC/m3/d
        # *layer_thickness*ndays -> to mgC/m2

        POC_resp_mgC_m2_vs_depth_list[:,iRespi] = PARR_list.copy()*list_dens_PARR*Oxy2C*mol2gC*layer_thickness*ndays/1000
        POC_resp_mgC_m2_std_vs_depth_list[:,iRespi] = PARR_std_list.copy()*list_dens_PARR*Oxy2C*mol2gC*layer_thickness*ndays/1000

    ########################################################################################################################
    # Here I calculate the carbon budget for depth0—depthf layer
    ########################################################################################################################
    Date_Num_Flux = x_filtered
    depth_POC_resp = list_depth_PARR
    depth_02_resp = depth_isopycnal

    ############### I calculate the integrated POC (MiP+MaP+bbp), between depth0 and depthf, for day0 and dayf. I transform it to mgC/m2

    # I extract the index of Integrated_POC_mgC_m3 which correspond to day0 (and dayf)
    tmp = list_dates_Integrated_POC - day0_float
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = list_dates_Integrated_POC - dayf_float
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]

    Integrated_POC_day0_mgC_m2 = Integrated_POC_mgC_m3[idx0] * layer_thickness
    Integrated_POC_dayf_mgC_m2 = Integrated_POC_mgC_m3[idxf] * layer_thickness
    Integrated_POC_day0_mgC_m2_std = Integrated_POC_mgC_m3_std[idx0] * layer_thickness
    Integrated_POC_dayf_mgC_m2_std = Integrated_POC_mgC_m3_std[idxf] * layer_thickness
    Delta_Integrated_POC = Integrated_POC_dayf_mgC_m2 - Integrated_POC_day0_mgC_m2
    Delta_Integrated_POC_std = np.sqrt( Integrated_POC_dayf_mgC_m2_std**2 + Integrated_POC_day0_mgC_m2_std**2 )

    Integrated_POC_extended_day0_mgC_m2 = Integrated_POC_extended_mgC_m3[idx0] * layer_thickness
    Integrated_POC_extended_dayf_mgC_m2 = Integrated_POC_extended_mgC_m3[idxf] * layer_thickness
    Integrated_POC_extended_day0_mgC_m2_std = Integrated_POC_extended_mgC_m3_std[idx0] * layer_thickness
    Integrated_POC_extended_dayf_mgC_m2_std = Integrated_POC_extended_mgC_m3_std[idxf] * layer_thickness
    Delta_Integrated_POC_extended = Integrated_POC_extended_dayf_mgC_m2 - Integrated_POC_extended_day0_mgC_m2
    Delta_Integrated_POC_extended_std = np.sqrt( Integrated_POC_extended_dayf_mgC_m2_std**2 + Integrated_POC_extended_day0_mgC_m2_std**2 )

    ############### I calculate the amount of POC entering from depht0 and exiting from dayf between day0 and dayf (in mgC/m2)

    # I extract the index of Flux_depth0/Flux_depthf which correspond to day0 (and dayf)
    tmp = Date_Num_Flux - day0_float
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = Date_Num_Flux - dayf_float
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]

    tmp = Flux_depth0[idx0:idxf]
    Flux_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_depthf[idx0:idxf]
    Flux_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_eta_b_depth0[idx0:idxf]
    Flux_eta_b_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_eta_b_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_eta_b_depthf[idx0:idxf]
    Flux_eta_b_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_eta_b_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_extended_depth0[idx0:idxf]
    Flux_extended_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_extended_depthf[idx0:idxf]
    Flux_extended_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_depthf_mgC_m2_std = np.std(tmp) * ndays

    tmp = Flux_extended_eta_b_depth0[idx0:idxf]
    Flux_extended_eta_b_depth0_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_eta_b_depth0_mgC_m2_std = np.std(tmp) * ndays
    tmp = Flux_extended_eta_b_depthf[idx0:idxf]
    Flux_extended_eta_b_depthf_mgC_m2 = np.mean(tmp) * ndays
    Flux_extended_eta_b_depthf_mgC_m2_std = np.std(tmp) * ndays

    Delta_flux = Flux_depth0_mgC_m2 - Flux_depthf_mgC_m2
    Delta_flux_eta_b = Flux_eta_b_depth0_mgC_m2 - Flux_eta_b_depthf_mgC_m2
    Delta_flux_extended = Flux_extended_depth0_mgC_m2 - Flux_extended_depthf_mgC_m2
    Delta_flux_extended_eta_b = Flux_extended_eta_b_depth0_mgC_m2 - Flux_extended_eta_b_depthf_mgC_m2

    Delta_flux_std = np.sqrt( Flux_depth0_mgC_m2_std**2 + Flux_depthf_mgC_m2_std**2 )
    Delta_flux_eta_b_std = np.sqrt( Flux_eta_b_depth0_mgC_m2_std**2 + Flux_eta_b_depthf_mgC_m2_std**2 )
    Delta_flux_extended_std = np.sqrt( Flux_extended_depth0_mgC_m2_std**2 + Flux_extended_depthf_mgC_m2_std**2 )
    Delta_flux_extended_eta_b_std = np.sqrt( Flux_extended_eta_b_depth0_mgC_m2_std**2 + Flux_extended_eta_b_depthf_mgC_m2_std**2 )

    Theoretical_Budget = Delta_flux - Delta_Integrated_POC
    Theoretical_Budget_eta_b = Delta_flux_eta_b - Delta_Integrated_POC
    Theoretical_Budget_extended = Delta_flux_extended - Delta_Integrated_POC_extended
    Theoretical_Budget_extended_eta_b = Delta_flux_extended_eta_b - Delta_Integrated_POC_extended

    Theoretical_Budget_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_eta_b_std = np.sqrt( Delta_flux_eta_b_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_extended_std**2 )
    Theoretical_Budget_extended_eta_b_std = np.sqrt( Delta_flux_extended_eta_b_std**2 + Delta_Integrated_POC_extended_std**2 )

    ############### I calculate the PARR and oxygen consumption between depth0 and depthf (in mgC/m2)

    # I extract the index of POC_resp_mgC_m2_vs_depth which correspond to depth0 (and depthf)
    tmp = depth_POC_resp - depth0
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = depth_POC_resp - depthf
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    POC_resp_mgC_m2_list=np.mean(POC_resp_mgC_m2_vs_depth_list[idx0:idxf,:],0)
    POC_resp_mgC_m2_std_list=np.mean(POC_resp_mgC_m2_std_vs_depth_list[idx0:idxf,:],0)

    # I extract the index of O2_resp_mgC_m2_vs_depth which correspond to depth0 (and depthf)
    tmp = depth_02_resp - depth0
    idx0 = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    tmp = depth_02_resp - depthf
    idxf = np.where(np.abs(tmp) == (np.abs(tmp)).min())[0][0]
    O2_resp_mgC_m2=np.mean(O2_resp_mgC_m2_vs_depth[idx0:idxf])*-1
    O2_resp_mgC_m2_ci=np.mean(O2_resp_mgC_m2_vs_depth_ci[idx0:idxf,:],0)*-1

    ############### I return the data
    return Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_eta_b,Theoretical_Budget_eta_b_std,\
           Theoretical_Budget_extended,Theoretical_Budget_extended_std,\
           Theoretical_Budget_extended_eta_b,Theoretical_Budget_extended_eta_b_std,\
           POC_resp_mgC_m2_list,POC_resp_mgC_m2_std_list,O2_resp_mgC_m2,O2_resp_mgC_m2_ci,list_Respi_types,n_profiles

########################################################################################################################
########################################################################################################################
########################################################################################################################
# End of function definition: I start the loop
########################################################################################################################
########################################################################################################################
########################################################################################################################

ndays_list=np.r_[ndays1:ndays2:delta_days];ndays_list=np.append(ndays_list,ndays2);ndays_list=np.unique(ndays_list)

ndays=ndays_list[0]
for ndays in ndays_list:
    dayf1 = day0 + timedelta(days=int(ndays))  # final date for the carbon budget calculation

    Theoretical_Budget_list_w1 = np.array([])
    Theoretical_Budget_eta_b_list_w1 = np.array([])
    Theoretical_Budget_extended_list_w1 = np.array([])
    Theoretical_Budget_extended_eta_b_list_w1 = np.array([])
    Theoretical_Budget_std_list_w1 = np.array([])
    Theoretical_Budget_eta_b_std_list_w1 = np.array([])
    Theoretical_Budget_extended_std_list_w1 = np.array([])
    Theoretical_Budget_extended_eta_b_std_list_w1 = np.array([])
    POC_resp_mgC_m2_list_w1 = np.array([])
    POC_resp_mgC_m2_std_list_w1 = np.array([])
    O2_resp_mgC_m2_list_w1 = np.array([])
    O2_resp_mgC_m2_ci_list_w1 = np.array([])
    Theoretical_Budget_list_w2 = np.array([])
    Theoretical_Budget_eta_b_list_w2 = np.array([])
    Theoretical_Budget_extended_list_w2 = np.array([])
    Theoretical_Budget_extended_eta_b_list_w2 = np.array([])
    Theoretical_Budget_std_list_w2 = np.array([])
    Theoretical_Budget_eta_b_std_list_w2 = np.array([])
    Theoretical_Budget_extended_std_list_w2 = np.array([])
    Theoretical_Budget_extended_eta_b_std_list_w2 = np.array([])
    POC_resp_mgC_m2_list_w2 = np.array([])
    POC_resp_mgC_m2_std_list_w2 = np.array([])
    O2_resp_mgC_m2_list_w2 = np.array([])
    O2_resp_mgC_m2_ci_list_w2 = np.array([])
    depth0=depth0_list[0]
    for depth0 in depth0_list:
        depthf = depth0+layer_thickness
        (Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_eta_b,Theoretical_Budget_eta_b_std,Theoretical_Budget_extended,Theoretical_Budget_extended_std,Theoretical_Budget_extended_eta_b,Theoretical_Budget_extended_eta_b_std,POC_resp_mgC_m2,POC_resp_mgC_m2_std,O2_resp_mgC_m2,O2_resp_mgC_m2_ci,RespirationTypes,n_profiles)=carbon_budget_calculation(depth0,depthf,day0,dayf1)
        Theoretical_Budget_list_w1=np.append(Theoretical_Budget_list_w1,Theoretical_Budget)
        Theoretical_Budget_eta_b_list_w1=np.append(Theoretical_Budget_eta_b_list_w1,Theoretical_Budget)
        Theoretical_Budget_extended_list_w1=np.append(Theoretical_Budget_extended_list_w1,Theoretical_Budget_extended)
        Theoretical_Budget_extended_eta_b_list_w1=np.append(Theoretical_Budget_extended_eta_b_list_w1,Theoretical_Budget_extended_eta_b)
        Theoretical_Budget_std_list_w1=np.append(Theoretical_Budget_std_list_w1,Theoretical_Budget)
        Theoretical_Budget_eta_b_std_list_w1=np.append(Theoretical_Budget_eta_b_std_list_w1,Theoretical_Budget_std)
        Theoretical_Budget_extended_std_list_w1=np.append(Theoretical_Budget_extended_std_list_w1,Theoretical_Budget_extended_std)
        Theoretical_Budget_extended_eta_b_std_list_w1=np.append(Theoretical_Budget_extended_eta_b_std_list_w1,Theoretical_Budget_extended_eta_b_std)
        POC_resp_mgC_m2_list_w1=np.append(POC_resp_mgC_m2_list_w1,POC_resp_mgC_m2,axis=0)
        POC_resp_mgC_m2_std_list_w1=np.append(POC_resp_mgC_m2_std_list_w1,POC_resp_mgC_m2_std,axis=0)
        O2_resp_mgC_m2_list_w1=np.append(O2_resp_mgC_m2_list_w1,O2_resp_mgC_m2)
        O2_resp_mgC_m2_ci_list_w1=np.append(O2_resp_mgC_m2_ci_list_w1,O2_resp_mgC_m2_ci,axis=0)
        (Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_eta_b,Theoretical_Budget_eta_b_std,Theoretical_Budget_extended,Theoretical_Budget_extended_std,Theoretical_Budget_extended_eta_b,Theoretical_Budget_extended_eta_b_std,POC_resp_mgC_m2,POC_resp_mgC_m2_std,O2_resp_mgC_m2,O2_resp_mgC_m2_ci,RespirationTypes,n_profiles)=carbon_budget_calculation(depth0,depthf,dayf1,dayf)
        Theoretical_Budget_list_w2=np.append(Theoretical_Budget_list_w2,Theoretical_Budget)
        Theoretical_Budget_eta_b_list_w2=np.append(Theoretical_Budget_eta_b_list_w2,Theoretical_Budget)
        Theoretical_Budget_extended_list_w2=np.append(Theoretical_Budget_extended_list_w2,Theoretical_Budget_extended)
        Theoretical_Budget_extended_eta_b_list_w2=np.append(Theoretical_Budget_extended_eta_b_list_w2,Theoretical_Budget_extended_eta_b)
        Theoretical_Budget_std_list_w2=np.append(Theoretical_Budget_std_list_w2,Theoretical_Budget)
        Theoretical_Budget_eta_b_std_list_w2=np.append(Theoretical_Budget_eta_b_std_list_w2,Theoretical_Budget_std)
        Theoretical_Budget_extended_std_list_w2=np.append(Theoretical_Budget_extended_std_list_w2,Theoretical_Budget_extended_std)
        Theoretical_Budget_extended_eta_b_std_list_w2=np.append(Theoretical_Budget_extended_eta_b_std_list_w2,Theoretical_Budget_extended_eta_b_std)
        POC_resp_mgC_m2_list_w2=np.append(POC_resp_mgC_m2_list_w2,POC_resp_mgC_m2,axis=0)
        POC_resp_mgC_m2_std_list_w2=np.append(POC_resp_mgC_m2_std_list_w2,POC_resp_mgC_m2_std,axis=0)
        O2_resp_mgC_m2_list_w2=np.append(O2_resp_mgC_m2_list_w2,O2_resp_mgC_m2)
        O2_resp_mgC_m2_ci_list_w2=np.append(O2_resp_mgC_m2_ci_list_w2,O2_resp_mgC_m2_ci,axis=0)

    O2_resp_mgC_m2_ci_list_w1=O2_resp_mgC_m2_ci_list_w1.reshape(depth0_list.size,2)
    POC_resp_mgC_m2_list_w1=POC_resp_mgC_m2_list_w1.reshape(depth0_list.size,len(RespirationTypes))
    POC_resp_mgC_m2_std_list_w1=POC_resp_mgC_m2_std_list_w1.reshape(depth0_list.size,len(RespirationTypes))
    O2_resp_mgC_m2_ci_list_w2=O2_resp_mgC_m2_ci_list_w2.reshape(depth0_list.size,2)
    POC_resp_mgC_m2_list_w2=POC_resp_mgC_m2_list_w2.reshape(depth0_list.size,len(RespirationTypes))
    POC_resp_mgC_m2_std_list_w2=POC_resp_mgC_m2_std_list_w2.reshape(depth0_list.size,len(RespirationTypes))

    ##################################################################
    # I plot
    fs=10
    width, height = 0.78, 0.75

    ##################################################################
    # First plot: budget calculated without considering smallest size classes and with new eta and b values, time window 1
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height])
    plt.plot(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, 'k')
    plt.scatter(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, c='black',s=5)
    plt.fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(2,3):
        plt.plot(POC_resp_mgC_m2_list_w1[:,iResp], depth0_list + layer_thickness / 2, c='b')

    plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,3], POC_resp_mgC_m2_list[:,4], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher)')

    plt.plot(POC_resp_mgC_m2_list_w1[:, 0], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 3], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 4], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 5], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 6], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    plt.plot(Theoretical_Budget_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red')
    plt.scatter(Theoretical_Budget_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red', s=5)
    plt.fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_eta_b_list_w1 - Theoretical_Budget_eta_b_std_list_w1*0.5, Theoretical_Budget_eta_b_list_w1 + Theoretical_Budget_eta_b_std_list_w1*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')

    plt.xlim(-200,7500)
    plt.ylabel('Depth (m)', fontsize=fs)
    plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.title('No small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (day0.year,day0.month,day0.day,dayf1.year,dayf1.month,dayf1.day), fontsize=9)
    plt.legend(fontsize=7)
    plt.gca().invert_yaxis()
    ax.text(-0.05, 1.125, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an34/CarbonBudget_vs_depth_IMday%d%02d%02d_TW1_from%d%02d%02d_0101eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,day0.year,day0.month,day0.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Second plot: budget calculated without considering smallest size classes and with new eta and b values, time window 2
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height])
    plt.plot(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, 'k')
    plt.scatter(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, c='black',s=5)
    plt.fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w2[:, 1], O2_resp_mgC_m2_ci_list_w2[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(2,3):
        plt.plot(POC_resp_mgC_m2_list_w2[:,iResp], depth0_list + layer_thickness / 2, c='b')

    plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w2[:,iResp]-POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w2[:,iResp]+POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,3], POC_resp_mgC_m2_list[:,4], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher)')

    plt.plot(POC_resp_mgC_m2_list_w2[:, 0], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 3], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 4], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 5], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 6], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    plt.plot(Theoretical_Budget_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red')
    plt.scatter(Theoretical_Budget_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red', s=5)
    plt.fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_eta_b_list_w2 - Theoretical_Budget_eta_b_std_list_w2*0.5, Theoretical_Budget_eta_b_list_w2 + Theoretical_Budget_eta_b_std_list_w2*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')

    plt.xlim(-200,7500)
    plt.ylabel('Depth (m)', fontsize=fs)
    plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.title('No small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day), fontsize=9)
    plt.legend(fontsize=7)
    plt.gca().invert_yaxis()
    ax.text(-0.05, 1.125, 'c', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an34/CarbonBudget_vs_depth_IMday%d%02d%02d_TW2_to%d%02d%02d_0101eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Third plot: budget calculated without considering smallest size classes and with new eta and b values, time window 1
    # I plot also the significant of the difference between the PARR, oxygen, and theo POC respired
    iResp=2
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w1[:,iResp]-Theoretical_Budget_eta_b_list_w1) / np.sqrt( POC_resp_mgC_m2_std_list_w1[:,iResp]**2 + Theoretical_Budget_eta_b_std_list_w1**2 )
    PARR_Theo = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w1[:,iResp]-O2_resp_mgC_m2_list_w1) / np.sqrt( POC_resp_mgC_m2_std_list_w1[:,iResp]**2 + (O2_resp_mgC_m2_ci_list_w1[:, 0]-O2_resp_mgC_m2_ci_list_w1[:, 1])**2 )
    PARR_Oxy = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(Theoretical_Budget_eta_b_list_w1-O2_resp_mgC_m2_list_w1) / np.sqrt( Theoretical_Budget_eta_b_std_list_w1 + (O2_resp_mgC_m2_ci_list_w1[:, 0]-O2_resp_mgC_m2_ci_list_w1[:, 1])**2 )
    Theo_Oxy = 2 * (norm.sf(twoSamplesZtest))

    fs=9
    width2 = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(3.5,3.5), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
    plt.subplots_adjust(wspace=0.05)
    ax[0].set_position([0.18, 0.15, width2, height])
    ax[1].set_position([0.18+width2+0.05, 0.15, 1-(0.18+width2+0.1), height])
    ax[0].plot(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, 'k')
    ax[0].scatter(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, c='black',s=5)
    ax[0].fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(2,3):
        ax[0].plot(POC_resp_mgC_m2_list_w1[:,iResp], depth0_list + layer_thickness / 2, c='b')

    ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,3], POC_resp_mgC_m2_list[:,4], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 0], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 3], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 4], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 5], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 6], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    ax[0].plot(Theoretical_Budget_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red')
    ax[0].scatter(Theoretical_Budget_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red', s=5)
    ax[0].fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_eta_b_list_w1 - Theoretical_Budget_eta_b_std_list_w1*0.5, Theoretical_Budget_eta_b_list_w1 + Theoretical_Budget_eta_b_std_list_w1*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')
    ax[0].set_xlim(-200,7500)
    ax[0].set_ylabel('Depth (m)', fontsize=fs)
    ax[0].set_xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.suptitle('No small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (day0.year,day0.month,day0.day,dayf1.year,dayf1.month,dayf1.day), fontsize=9)
    ax[0].legend(fontsize=7)
    ax[0].invert_yaxis()
    ax[0].text(-0.2, 1.115, 'b', transform=ax[0].transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
    filled_marker_style = dict(marker='o',color='white',linewidth=0.00001, markersize=6, markeredgecolor='black',fillstyle='left',markeredgewidth=.3)
    ax[1].plot(1*np.ones((sum(PARR_Theo>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Theo>0.05],markerfacecolor='blue', markerfacecoloralt='red',**filled_marker_style)
    ax[1].plot(2*np.ones((sum(PARR_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Oxy>0.05],markerfacecolor='blue', markerfacecoloralt='gray',**filled_marker_style)
    ax[1].plot(3*np.ones((sum(Theo_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[Theo_Oxy>0.05],markerfacecolor='gray', markerfacecoloralt='red',**filled_marker_style)
    ax[1].set_xticks(np.array([0.5, 1.5, 2.5, 3.5]))
    ax[1].set_xticklabels([])
    ax[1].grid(color='k', linestyle='dashed', linewidth=0.5, which='both')
    ax[1].text(-0.2, -0.040, 'Statistical\ndifference', transform=ax[1].transAxes, fontsize=9, va='top',
               ha='left')  # ,fontfamily='helvetica'
    plt.savefig('../Plots/an34/WithStatisticalDifference/CarbonBudget_vs_depth_IMday%d%02d%02d_TW1_from%d%02d%02d_0201eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,day0.year,day0.month,day0.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Fourth plot: budget calculated without considering smallest size classes and with new eta and b values, time window 2
    # I plot also the significant of the difference between the PARR, oxygen, and theo POC respired
    iResp=2
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w2[:,iResp]-Theoretical_Budget_eta_b_list_w2) / np.sqrt( POC_resp_mgC_m2_std_list_w2[:,iResp]**2 + Theoretical_Budget_eta_b_std_list_w2**2 )
    PARR_Theo = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w2[:,iResp]-O2_resp_mgC_m2_list_w2) / np.sqrt( POC_resp_mgC_m2_std_list_w2[:,iResp]**2 + (O2_resp_mgC_m2_ci_list_w2[:, 0]-O2_resp_mgC_m2_ci_list_w2[:, 1])**2 )
    PARR_Oxy = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(Theoretical_Budget_eta_b_list_w2-O2_resp_mgC_m2_list_w2) / np.sqrt( Theoretical_Budget_eta_b_std_list_w2 + (O2_resp_mgC_m2_ci_list_w2[:, 0]-O2_resp_mgC_m2_ci_list_w2[:, 1])**2 )
    Theo_Oxy = 2 * (norm.sf(twoSamplesZtest))

    fs=9
    width2 = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(3.5,3.5), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
    plt.subplots_adjust(wspace=0.05)
    ax[0].set_position([0.18, 0.15, width2, height])
    ax[1].set_position([0.18+width2+0.05, 0.15, 1-(0.18+width2+0.1), height])
    ax[0].plot(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, 'k')
    ax[0].scatter(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, c='black',s=5)
    ax[0].fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w2[:, 1], O2_resp_mgC_m2_ci_list_w2[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(2,3):
        ax[0].plot(POC_resp_mgC_m2_list_w2[:,iResp], depth0_list + layer_thickness / 2, c='b')

    ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w2[:,iResp]-POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w2[:,iResp]+POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,3], POC_resp_mgC_m2_list[:,4], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 0], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 3], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 4], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 5], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 6], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    ax[0].plot(Theoretical_Budget_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red')
    ax[0].scatter(Theoretical_Budget_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red', s=5)
    ax[0].fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_eta_b_list_w2 - Theoretical_Budget_eta_b_std_list_w2*0.5, Theoretical_Budget_eta_b_list_w2 + Theoretical_Budget_eta_b_std_list_w2*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')
    ax[0].set_xlim(-200,7500)
    ax[0].set_ylabel('Depth (m)', fontsize=fs)
    ax[0].set_xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.suptitle('No small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day), fontsize=9)
    ax[0].legend(fontsize=7)
    ax[0].invert_yaxis()
    ax[0].text(-0.2, 1.115, 'b', transform=ax[0].transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
    filled_marker_style = dict(marker='o',color='white',linewidth=0.00001, markersize=6, markeredgecolor='black',fillstyle='left',markeredgewidth=.3)
    ax[1].plot(1*np.ones((sum(PARR_Theo>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Theo>0.05],markerfacecolor='blue', markerfacecoloralt='red',**filled_marker_style)
    ax[1].plot(2*np.ones((sum(PARR_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Oxy>0.05],markerfacecolor='blue', markerfacecoloralt='gray',**filled_marker_style)
    ax[1].plot(3*np.ones((sum(Theo_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[Theo_Oxy>0.05],markerfacecolor='gray', markerfacecoloralt='red',**filled_marker_style)
    ax[1].set_xticks(np.array([0.5, 1.5, 2.5, 3.5]))
    ax[1].set_xticklabels([])
    ax[1].grid(color='k', linestyle='dashed', linewidth=0.5, which='both')
    ax[1].text(-0.2, -0.040, 'Statistical\ndifference', transform=ax[1].transAxes, fontsize=9, va='top',
               ha='left')  # ,fontfamily='helvetica'
    plt.savefig('../Plots/an34/WithStatisticalDifference/CarbonBudget_vs_depth_IMday%d%02d%02d_TW2_to%d%02d%02d_0201eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Fifth plot: budget calculated considering smallest size classes and with different eta and b values (based on the
    # relationship between sinking speed and particle size), time window 1
    fs=10
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height])
    plt.plot(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, 'k')
    plt.scatter(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, c='black',s=5)
    plt.fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(9,10):
        plt.plot(POC_resp_mgC_m2_list_w1[:,iResp], depth0_list + layer_thickness / 2, c='b')

    plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,10], POC_resp_mgC_m2_list[:,11], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher)')

    plt.plot(POC_resp_mgC_m2_list_w1[:, 7], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 10], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 11], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 12], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    plt.plot(POC_resp_mgC_m2_list_w1[:, 13], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    plt.plot(Theoretical_Budget_extended_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red')
    plt.scatter(Theoretical_Budget_extended_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red', s=5)
    plt.fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_extended_eta_b_list_w1 - Theoretical_Budget_extended_eta_b_std_list_w1*0.5, Theoretical_Budget_extended_eta_b_list_w1 + Theoretical_Budget_extended_eta_b_std_list_w1*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')


    plt.xlim(-200,7500)
    plt.ylabel('Depth (m)', fontsize=fs)
    plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.title('With small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (day0.year,day0.month,day0.day,dayf1.year,dayf1.month,dayf1.day), fontsize=9)
    plt.legend(fontsize=7)
    plt.gca().invert_yaxis()
    ax.text(-0.05, 1.125, 'd', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an34/CarbonBudget_vs_depth_IMday%d%02d%02d_TW1_from%d%02d%02d_0102extended_eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,day0.year,day0.month,day0.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Sixth plot: budget calculated considering smallest size classes and with different eta and b values (based on the
    # relationship between sinking speed and particle size), time window 2
    fs=10
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.18, 0.15, width, height])
    plt.plot(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, 'k')
    plt.scatter(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, c='black',s=5)
    plt.fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w2[:, 1], O2_resp_mgC_m2_ci_list_w2[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(9,10):
        plt.plot(POC_resp_mgC_m2_list_w2[:,iResp], depth0_list + layer_thickness / 2, c='b')

    plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w2[:,iResp]-POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w2[:,iResp]+POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # plt.fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,10], POC_resp_mgC_m2_list[:,11], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher)')

    plt.plot(POC_resp_mgC_m2_list_w2[:, 7], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 10], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 11], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 12], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    plt.plot(POC_resp_mgC_m2_list_w2[:, 13], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    plt.plot(Theoretical_Budget_extended_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red')
    plt.scatter(Theoretical_Budget_extended_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red', s=5)
    plt.fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_extended_eta_b_list_w2 - Theoretical_Budget_extended_eta_b_std_list_w2*0.5, Theoretical_Budget_extended_eta_b_list_w2 + Theoretical_Budget_extended_eta_b_std_list_w2*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')


    plt.xlim(-200,7500)
    plt.ylabel('Depth (m)', fontsize=fs)
    plt.xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.title('With small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day), fontsize=9)
    plt.legend(fontsize=7)
    plt.gca().invert_yaxis()
    ax.text(-0.05, 1.125, 'e', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an34/CarbonBudget_vs_depth_IMday%d%02d%02d_TW2_to%d%02d%02d_0102extended_eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Seventh plot: budget calculated considering smallest size classes and with different eta and b values (based on the
    # relationship between sinking speed and particle size), time window 1
    # I plot also the significant of the difference between the PARR, oxygen, and theo POC respired
    iResp=9
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w1[:,iResp]-Theoretical_Budget_extended_eta_b_list_w1) / np.sqrt( POC_resp_mgC_m2_std_list_w1[:,iResp]**2 + Theoretical_Budget_extended_eta_b_std_list_w1**2 )
    PARR_Theo = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w1[:,iResp]-O2_resp_mgC_m2_list_w1) / np.sqrt( POC_resp_mgC_m2_std_list_w1[:,iResp]**2 + (O2_resp_mgC_m2_ci_list_w1[:, 0]-O2_resp_mgC_m2_ci_list_w1[:, 1])**2 )
    PARR_Oxy = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(Theoretical_Budget_extended_eta_b_list_w1-O2_resp_mgC_m2_list_w1) / np.sqrt( Theoretical_Budget_extended_eta_b_std_list_w1 + (O2_resp_mgC_m2_ci_list_w1[:, 0]-O2_resp_mgC_m2_ci_list_w1[:, 1])**2 )
    Theo_Oxy = 2 * (norm.sf(twoSamplesZtest))

    fs=9
    width2 = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(3.5,3.5), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
    plt.subplots_adjust(wspace=0.05)
    ax[0].set_position([0.18, 0.15, width2, height])
    ax[1].set_position([0.18+width2+0.05, 0.15, 1-(0.18+width2+0.1), height])
    ax[0].plot(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, 'k')
    ax[0].scatter(O2_resp_mgC_m2_list_w1,depth0_list+layer_thickness/2, c='black',s=5)
    ax[0].fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w1[:, 1], O2_resp_mgC_m2_ci_list_w1[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(9,10):
        ax[0].plot(POC_resp_mgC_m2_list_w1[:,iResp], depth0_list + layer_thickness / 2, c='b')

    ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w1[:,iResp]-POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w1[:,iResp]+POC_resp_mgC_m2_std_list_w1[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,10], POC_resp_mgC_m2_list[:,11], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 7], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 10], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 11], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 12], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    ax[0].plot(POC_resp_mgC_m2_list_w1[:, 13], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    ax[0].plot(Theoretical_Budget_extended_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red')
    ax[0].scatter(Theoretical_Budget_extended_eta_b_list_w1, depth0_list + layer_thickness / 2, c='red', s=5)
    ax[0].fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_extended_eta_b_list_w1 - Theoretical_Budget_extended_eta_b_std_list_w1*0.5, Theoretical_Budget_extended_eta_b_list_w1 + Theoretical_Budget_extended_eta_b_std_list_w1*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')
    ax[0].set_xlim(-200,7500)
    ax[0].set_ylabel('Depth (m)', fontsize=fs)
    ax[0].set_xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.suptitle('With small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (day0.year,day0.month,day0.day,dayf1.year,dayf1.month,dayf1.day), fontsize=9)
    ax[0].legend(fontsize=7)
    ax[0].invert_yaxis()
    ax[0].text(-0.2, 1.115, 'b', transform=ax[0].transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
    filled_marker_style = dict(marker='o',color='white',linewidth=0.00001, markersize=6, markeredgecolor='black',fillstyle='left',markeredgewidth=.3)
    ax[1].plot(1*np.ones((sum(PARR_Theo>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Theo>0.05],markerfacecolor='blue', markerfacecoloralt='red',**filled_marker_style)
    ax[1].plot(2*np.ones((sum(PARR_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Oxy>0.05],markerfacecolor='blue', markerfacecoloralt='gray',**filled_marker_style)
    ax[1].plot(3*np.ones((sum(Theo_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[Theo_Oxy>0.05],markerfacecolor='gray', markerfacecoloralt='red',**filled_marker_style)
    ax[1].set_xticks(np.array([0.5, 1.5, 2.5, 3.5]))
    ax[1].set_xticklabels([])
    ax[1].grid(color='k', linestyle='dashed', linewidth=0.5, which='both')
    ax[1].text(-0.2, -0.040, 'Statistical\ndifference', transform=ax[1].transAxes, fontsize=9, va='top',
               ha='left')  # ,fontfamily='helvetica'
    plt.savefig('../Plots/an34/WithStatisticalDifference/CarbonBudget_vs_depth_IMday%d%02d%02d_TW1_from%d%02d%02d_0202extended_eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,day0.year,day0.month,day0.day) ,dpi=200)
    plt.close()

    ##################################################################
    # Eighth plot: budget calculated considering smallest size classes and with different eta and b values (based on the
    # relationship between sinking speed and particle size), time window 2
    # I plot also the significant of the difference between the PARR, oxygen, and theo POC respired
    iResp=9
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w2[:,iResp]-Theoretical_Budget_extended_eta_b_list_w2) / np.sqrt( POC_resp_mgC_m2_std_list_w2[:,iResp]**2 + Theoretical_Budget_extended_eta_b_std_list_w2**2 )
    PARR_Theo = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(POC_resp_mgC_m2_list_w2[:,iResp]-O2_resp_mgC_m2_list_w2) / np.sqrt( POC_resp_mgC_m2_std_list_w2[:,iResp]**2 + (O2_resp_mgC_m2_ci_list_w2[:, 0]-O2_resp_mgC_m2_ci_list_w2[:, 1])**2 )
    PARR_Oxy = 2 * (norm.sf(twoSamplesZtest))
    twoSamplesZtest=np.sqrt(n_profiles)*abs(Theoretical_Budget_extended_eta_b_list_w2-O2_resp_mgC_m2_list_w2) / np.sqrt( Theoretical_Budget_extended_eta_b_std_list_w2 + (O2_resp_mgC_m2_ci_list_w2[:, 0]-O2_resp_mgC_m2_ci_list_w2[:, 1])**2 )
    Theo_Oxy = 2 * (norm.sf(twoSamplesZtest))

    fs=9
    width2 = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(3.5,3.5), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
    plt.subplots_adjust(wspace=0.05)
    ax[0].set_position([0.18, 0.15, width2, height])
    ax[1].set_position([0.18+width2+0.05, 0.15, 1-(0.18+width2+0.1), height])
    ax[0].plot(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, 'k')
    ax[0].scatter(O2_resp_mgC_m2_list_w2,depth0_list+layer_thickness/2, c='black',s=5)
    ax[0].fill_betweenx(depth0_list+layer_thickness/2, O2_resp_mgC_m2_ci_list_w2[:, 1], O2_resp_mgC_m2_ci_list_w2[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$')
    for iResp in range(9,10):
        ax[0].plot(POC_resp_mgC_m2_list_w2[:,iResp], depth0_list + layer_thickness / 2, c='b')

    ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list_w2[:,iResp]-POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5,
                      POC_resp_mgC_m2_list_w2[:,iResp]+POC_resp_mgC_m2_std_list_w2[:,iResp]*0.5, facecolor='b',
                      color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    # ax[0].fill_betweenx(depth0_list+layer_thickness/2, POC_resp_mgC_m2_list[:,10], POC_resp_mgC_m2_list[:,11], facecolor='b',
    #                     color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013;\nBelcher et al.)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 7], depth0_list + layer_thickness / 2, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 10], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 11], depth0_list + layer_thickness / 2, c='b',linestyle='dashed')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 12], depth0_list + layer_thickness / 2, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1)')
    ax[0].plot(POC_resp_mgC_m2_list_w2[:, 13], depth0_list + layer_thickness / 2, c='g',ls='-.',label='PARR\n($k_{rem}$=0.5)')
    ax[0].plot(Theoretical_Budget_extended_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red')
    ax[0].scatter(Theoretical_Budget_extended_eta_b_list_w2, depth0_list + layer_thickness / 2, c='red', s=5)
    ax[0].fill_betweenx(depth0_list + layer_thickness / 2, Theoretical_Budget_extended_eta_b_list_w2 - Theoretical_Budget_extended_eta_b_std_list_w2*0.5, Theoretical_Budget_extended_eta_b_list_w2 + Theoretical_Budget_extended_eta_b_std_list_w2*0.5,
                      facecolor='r', color='r', alpha=0.5, label='Theoretical\nPOC\nrespired')
    ax[0].set_xlim(-200,7500)
    ax[0].set_ylabel('Depth (m)', fontsize=fs)
    ax[0].set_xlabel('Carbon Consumption Rate (mgC/m$^2$)', fontsize=fs)
    plt.suptitle('With small size classes, eta=0.62,b=66\nStart date: %d-%02d-%02d; End date: %d-%02d-%02d' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day), fontsize=9)
    ax[0].legend(fontsize=7)
    ax[0].invert_yaxis()
    ax[0].text(-0.2, 1.115, 'b', transform=ax[0].transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
    ax[0].grid(color='k', linestyle='dashed', linewidth=0.5)
    filled_marker_style = dict(marker='o',color='white',linewidth=0.00001, markersize=6, markeredgecolor='black',fillstyle='left',markeredgewidth=.3)
    ax[1].plot(1*np.ones((sum(PARR_Theo>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Theo>0.05],markerfacecolor='blue', markerfacecoloralt='red',**filled_marker_style)
    ax[1].plot(2*np.ones((sum(PARR_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[PARR_Oxy>0.05],markerfacecolor='blue', markerfacecoloralt='gray',**filled_marker_style)
    ax[1].plot(3*np.ones((sum(Theo_Oxy>0.05),1)),(depth0_list+layer_thickness/2)[Theo_Oxy>0.05],markerfacecolor='gray', markerfacecoloralt='red',**filled_marker_style)
    ax[1].set_xticks(np.array([0.5, 1.5, 2.5, 3.5]))
    ax[1].set_xticklabels([])
    ax[1].grid(color='k', linestyle='dashed', linewidth=0.5, which='both')
    ax[1].text(-0.2, -0.040, 'Statistical\ndifference', transform=ax[1].transAxes, fontsize=9, va='top',
               ha='left')  # ,fontfamily='helvetica'
    plt.savefig('../Plots/an34/WithStatisticalDifference/CarbonBudget_vs_depth_IMday%d%02d%02d_TW2_to%d%02d%02d_0202extended_eta_b_an34.pdf' % (dayf1.year,dayf1.month,dayf1.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
    plt.close()



