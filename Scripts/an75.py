q10=3.4

########################################################################################################################
########################################################################################################################
########################################################################################################################
######### SUPPLEMENTARY FIG PARR BULK and OXY CONSUMPTION with different krem
########################################################################################################################
########################################################################################################################
########################################################################################################################
# region Supplementary Fig. PARR BULK and OXY CONSUMPTION with different krem
import numpy as np
import pandas as pd
import os,sys
import netCDF4 as nc
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
filename_coriolis='6903095_Sprof_all.nc'
########
import datetime,calendar
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from paruvpy import bbpPOC_to_bbp_abundance
from paruvpy import calculate_remin_func_abun
from paruvpy import ESD_limits_to_ESD_middle
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
import seawater as sw
import gsw
from lin_fit import lin_fit

#######################################################################
# I define the function to extract the POC content at the beginning and end of the time series, used in the function below
#######################################################################
def POC_day0dayf(t,POC):
    ndays=t.max()
    (interpol, slpe_ci, intrcp_ci, _, _) = lin_fit(t, POC)
    POC_day0=interpol.intercept / ndays# * layer_thickness
    POC_dayf=(interpol.slope*ndays+interpol.intercept) / ndays# * layer_thickness
    slpe_std=abs((slpe_ci - interpol.slope)[0]);intrcp_std=abs((intrcp_ci - interpol.intercept)[0])
    POC_day0_std = intrcp_std / ndays# * layer_thickness
    POC_dayf_std = np.sqrt(ndays**2*slpe_std**2+intrcp_std**2) / ndays# * layer_thickness
    return POC_day0,POC_dayf,POC_day0_std,POC_dayf_std

#######################################################################
# I define the function for the carbon budget calculation
#######################################################################
# region carbon_budget_calculation(dens0,densf,day0,dayf):
def carbon_budget_calculation(dens0,densf,day0,dayf):
    ########################################################################################################################
    # Starting parameters
    ########################################################################################################################
    # dayf = day0+timedelta(days=ndays) # final date for the carbon budget calculation
    ndays = (dayf - day0).days  # number of days
    delta_dens_flux = 0.025     # around of the density which I consider when extracting the flux
    Oxy2C = 0.89                # to convert from mol of oxygen to mol of carbon
    Oxy2C_std = Oxy2C * 0.4
    mol2gC = 12.0107            # to convert from mol of carbon to grams of carbon
    day0_float = calendar.timegm(day0.timetuple())
    dayf_float = calendar.timegm(dayf.timetuple())
    day0_datenum = matlab_datenum(day0.year,day0.month,day0.day,day0.hour,day0.minute,day0.second)
    dayf_datenum = matlab_datenum(dayf.year,dayf.month,dayf.day,dayf.hour,dayf.minute,dayf.second)

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
    dens=np.array(data['Potential density [kg/m3]'][sel_filename])
    Flux=np.array(data['Flux_mgC_m2'][sel_filename])
    # Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
    Flux_extended=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
    # Flux_extended_eta_b=np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
    MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
    MiP_POC_extended=np.array(data['Mip_POC_cont_mgC_m3_extendendTo0.0254sizeclass'][sel_filename])
    MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num=np.r_[0:Flux.size]
    for i in Date_Num:
        date_time_obj = datetime.datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    # I select the data only in the prescribed period
    list_dates=np.sort(np.unique(Date_Num))
    list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
    n_profiles=list_dates.size

    # I load the bbp data and I select only those in the prescribed period
    storedir = '%s/GIT/AC_Agulhas_eddy_2021/Data' % home
    a_file = open("%s/an18/data_an18.pkl" % storedir, "rb")
    data_an18 = pickle.load(a_file)
    bbp_POC = data_an18['bbp_POC']
    bbp_POC_Koestner = data_an18['bbp_POC_Koestner']
    Date_Num_bbp = data_an18['Date_Num_bbp']
    Date_Vec_bbp = data_an18['Date_Vec_bbp']
    depth_bbp = data_an18['depth_bbp']
    dens_bbp = data_an18['dens_bbp']
    temp_bbp = data_an18['temperature']
    a_file.close()

    sel_dates = (Date_Num_bbp>=day0_datenum)&(Date_Num_bbp<=dayf_datenum)
    Date_Num_bbp = Date_Num_bbp[sel_dates]
    Date_Vec_bbp = Date_Vec_bbp[sel_dates,:]
    depth_bbp = depth_bbp[sel_dates,:]
    dens_bbp = dens_bbp[sel_dates,:]
    temp_bbp = temp_bbp[sel_dates,:]
    bbp_POC = bbp_POC[sel_dates, :]
    bbp_POC_Koestner = bbp_POC_Koestner[sel_dates, :]
    # I calculate the bbp abundance
    sel=(bbp_POC==99999)&(bbp_POC_Koestner==99999)
    bbp_abun_tmp=bbpPOC_to_bbp_abundance(bbp_POC[~sel])
    bbp_abun=bbp_POC.copy()*0+99999;bbp_abun[~sel]=bbp_abun_tmp
    bbp_abun_Koestner_tmp=bbpPOC_to_bbp_abundance(bbp_POC_Koestner[~sel])
    bbp_abun_Koestner=bbp_POC.copy()*0+99999;bbp_abun_Koestner[~sel]=bbp_abun_Koestner_tmp
    # I calculate the bbp PARR in nmol_l_h
    esd_bbp = ESD_limits_to_ESD_middle(0.001, 0.025)
    bbp_PARR_tmp = calculate_remin_func_abun(bbp_abun[~sel], esd_bbp, temp_bbp[~sel], kRemPoc=0.013)
    bbp_PARR=bbp_POC.copy()*0+99999;bbp_PARR[~sel]=bbp_PARR_tmp
    bbp_PARR_Koestner_tmp = calculate_remin_func_abun(bbp_abun_Koestner[~sel], esd_bbp, temp_bbp[~sel], kRemPoc=0.013)
    bbp_PARR_Koestner=bbp_POC.copy()*0+99999;bbp_PARR_Koestner[~sel]=bbp_PARR_Koestner_tmp
    del bbp_abun_tmp,bbp_abun_Koestner_tmp,bbp_abun,bbp_abun_Koestner,bbp_PARR_tmp,bbp_PARR_Koestner_tmp
    # I convert the bpp PARR from nmol_l_h to micromol/kg/day
    bbp_PARR[~sel] = bbp_PARR[~sel] / 1000 * 24 / (dens_bbp[~sel] / 1000)
    bbp_PARR_Koestner[~sel] = bbp_PARR_Koestner[~sel] / 1000 * 24 / (dens_bbp[~sel] / 1000)

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_bbp_calendar = Date_Num_bbp.copy()
    for i in range(0,Date_Num_bbp_calendar.size):
        date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                                          Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
        Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
        # datetime.utcfromtimestamp(Date_Num[i])

    ########################################################################################################################
    # Here I calculate the integrated POC (i.e., MiP+MaP+bbp). To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because
    # some profiles only reach 400 m; (iii) I extract the mean MiP+MaP+bbp values between dens0 and densf and between day0 and
    # dayf (I obtain a time series); (iv) I calculate the bbp PARR between dens0 and densf and between day0 and dayf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
    MiP_extended_filtered=np.array([]);dens_MiP_extended_filtered=np.array([]);Date_Num_MiP_extended_filtered=np.array([])
    MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
    bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])
    bbp_Koestner_filtered=np.array([]);dens_bbp_Koestner_filtered=np.array([]);Date_Num_bbp_Koestner_filtered=np.array([])
    bbp_PARR_filtered=np.array([]);dens_bbp_PARR_filtered=np.array([]);Date_Num_bbp_PARR_filtered=np.array([])
    bbp_PARR_Koestner_filtered=np.array([]);dens_bbp_PARR_Koestner_filtered=np.array([]);Date_Num_bbp_PARR_Koestner_filtered=np.array([])

    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=MiP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_filtered = np.concatenate((MiP_filtered, z))
            Date_Num_MiP_filtered = np.concatenate((Date_Num_MiP_filtered, x2))
            dens_MiP_filtered = np.concatenate((dens_MiP_filtered, y2))
        z=MiP_POC_extended[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MiP_extended_filtered = np.concatenate((MiP_extended_filtered, z))
            Date_Num_MiP_extended_filtered = np.concatenate((Date_Num_MiP_extended_filtered, x2))
            dens_MiP_extended_filtered = np.concatenate((dens_MiP_extended_filtered, y2))
        z=MaP_POC[sel];x=Date_Num[sel];y=dens[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            MaP_filtered = np.concatenate((MaP_filtered, z))
            Date_Num_MaP_filtered = np.concatenate((Date_Num_MaP_filtered, x2))
            dens_MaP_filtered = np.concatenate((dens_MaP_filtered, y2))

    i=0
    for i in range(0, bbp_POC.shape[0]):
        #Cetinic
        z=bbp_POC[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>100] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_filtered = np.concatenate((bbp_filtered, z))
            Date_Num_bbp_filtered = np.concatenate((Date_Num_bbp_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_filtered = np.concatenate((dens_bbp_filtered, y2))
        #Cetinic PARR
        z=bbp_PARR[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>100] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_PARR_filtered = np.concatenate((bbp_PARR_filtered, z))
            Date_Num_bbp_PARR_filtered = np.concatenate((Date_Num_bbp_PARR_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_PARR_filtered = np.concatenate((dens_bbp_PARR_filtered, y2))
        #Koestner
        z=bbp_POC_Koestner[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>400] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_Koestner_filtered = np.concatenate((bbp_Koestner_filtered, z))
            Date_Num_bbp_Koestner_filtered = np.concatenate((Date_Num_bbp_Koestner_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_Koestner_filtered = np.concatenate((dens_bbp_Koestner_filtered, y2))
        #Koestner PARR
        z=bbp_PARR_Koestner[i,:];y=dens_bbp[i,:];x = Date_Num_bbp_calendar[i]
        z[z>400] = 99999
        sel2=(~np.isnan(z)) & (z != 99999);z=z[sel2];y2=y[sel2]
        sel3=z==0
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            z[sel3]=0
            bbp_PARR_Koestner_filtered = np.concatenate((bbp_PARR_Koestner_filtered, z))
            Date_Num_bbp_PARR_Koestner_filtered = np.concatenate((Date_Num_bbp_PARR_Koestner_filtered, np.tile(x,sum(sel2)) ))
            dens_bbp_PARR_Koestner_filtered = np.concatenate((dens_bbp_PARR_Koestner_filtered, y2))

    # I define the x and y arrays for the MiP+MaP+bbp interpolation
    x_filtered = np.linspace(Date_Num_bbp_filtered.min(), Date_Num_bbp_filtered.max(), ndays)
    y_filtered = np.linspace(dens_bbp_filtered.min(), dens_MaP_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    MiP_interp = griddata((Date_Num_MiP_filtered, dens_MiP_filtered), MiP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MiP_extended_interp = griddata((Date_Num_MiP_extended_filtered, dens_MiP_extended_filtered), MiP_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    MaP_interp = griddata((Date_Num_MaP_filtered, dens_MaP_filtered), MaP_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_interp = griddata((Date_Num_bbp_filtered, dens_bbp_filtered), bbp_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_Koestner_interp = griddata((Date_Num_bbp_Koestner_filtered, dens_bbp_Koestner_filtered), bbp_Koestner_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_PARR_interp = griddata((Date_Num_bbp_PARR_filtered, dens_bbp_PARR_filtered), bbp_PARR_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    bbp_PARR_Koestner_interp = griddata((Date_Num_bbp_PARR_Koestner_filtered, dens_bbp_PARR_Koestner_filtered), bbp_PARR_Koestner_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, I calculate the mean MiP+MaP+bbp (and std) between dens0 and densf between day0 and dayf
    sel_dens0_densf = (np.abs(y_filtered) >= dens0) & (np.abs(y_filtered) < densf)
    MiP_POC_dens0_densf=np.mean(MiP_interp[sel_dens0_densf,:],0)
    MiP_POC_extended_dens0_densf=np.mean(MiP_extended_interp[sel_dens0_densf,:],0)
    MaP_POC_dens0_densf=np.mean(MaP_interp[sel_dens0_densf,:],0)
    bbp_POC_dens0_densf=np.mean(bbp_interp[sel_dens0_densf,:],0)
    bbp_POC_Koestner_dens0_densf=np.mean(bbp_Koestner_interp[sel_dens0_densf,:],0)
    bbp_PARR_dens0_densf=np.mean(bbp_PARR_interp[sel_dens0_densf,:])
    bbp_PARR_Koestner_dens0_densf=np.mean(bbp_PARR_Koestner_interp[sel_dens0_densf,:])

    MiP_POC_dens0_densf_std = np.std(MiP_interp[sel_dens0_densf, :], 0)
    MiP_POC_extended_dens0_densf_std = np.std(MiP_extended_interp[sel_dens0_densf, :], 0)
    MaP_POC_dens0_densf_std = np.std(MaP_interp[sel_dens0_densf, :], 0)
    bbp_POC_dens0_densf_std = np.std(bbp_interp[sel_dens0_densf, :], 0)
    bbp_POC_Koestner_dens0_densf_std = np.std(bbp_Koestner_interp[sel_dens0_densf, :], 0)
    bbp_PARR_dens0_densf_std = np.std(bbp_PARR_interp[sel_dens0_densf, :])
    bbp_PARR_Koestner_dens0_densf_std = np.std(bbp_PARR_Koestner_interp[sel_dens0_densf, :])

    Integrated_POC_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf + bbp_POC_dens0_densf
    Integrated_POC_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    Integrated_POC_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_dens0_densf_std**2 )
    Integrated_POC_noBBP_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf
    Integrated_POC_noBBP_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf
    Integrated_POC_noBBP_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2  )
    Integrated_POC_noBBP_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 )
    Integrated_POC_Koestner_mgC_m3 = MiP_POC_dens0_densf + MaP_POC_dens0_densf + bbp_POC_Koestner_dens0_densf
    Integrated_POC_Koestner_extended_mgC_m3 = MiP_POC_extended_dens0_densf + MaP_POC_dens0_densf + bbp_POC_Koestner_dens0_densf
    Integrated_POC_Koestner_mgC_m3_std = np.sqrt( MiP_POC_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_Koestner_dens0_densf_std**2 )
    Integrated_POC_Koestner_extended_mgC_m3_std = np.sqrt( MiP_POC_extended_dens0_densf_std**2 + MaP_POC_dens0_densf_std**2 + bbp_POC_Koestner_dens0_densf_std**2 )

    ########################################################################################################################
    # Here I extract the flux values at dens0 and densf. To do so, (i) I filter it with a savgol function, then (ii) I
    # interpolate it over a regular grid in time and density. This step is necessary to have the flux at 600 m, because some
    # profiles only reach 400 m; (iii) I extract the flux values at dens0 and densf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    Flux_filtered=np.array([]);dens_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
    Flux_extended_filtered=np.array([]);dens_Flux_extended_filtered=np.array([]);Date_Num_Flux_extended_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=Flux[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_filtered = np.concatenate((Flux_filtered, z))
            Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
            dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y2))
        z=Flux_extended[sel];x=Date_Num[sel];y=dens[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_extended_filtered = np.concatenate((Flux_extended_filtered, z))
            Date_Num_Flux_extended_filtered = np.concatenate((Date_Num_Flux_extended_filtered, x2))
            dens_Flux_extended_filtered = np.concatenate((dens_Flux_extended_filtered, y2))

    # I define the x and y arrays for the Flux interpolation
    x_filtered = np.linspace(Date_Num_Flux_filtered.min(), Date_Num_Flux_filtered.max(), 100)
    y_filtered = np.linspace(dens_Flux_filtered.min(), dens_Flux_filtered.max(), 1000)
    x_filtered_g, y_filtered_g = np.meshgrid(x_filtered, y_filtered)
    # I interpolate
    Flux_interp = griddata((Date_Num_Flux_filtered, dens_Flux_filtered), Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")
    Flux_extended_interp = griddata((Date_Num_Flux_extended_filtered, dens_Flux_extended_filtered), Flux_extended_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, flux extraction at dens0 and densf

    sel_layer = (np.abs(y_filtered) >= dens0-delta_dens_flux) & (np.abs(y_filtered) < dens0+delta_dens_flux)
    Flux_dens0 = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_extended_dens0 = np.mean(Flux_extended_interp[sel_layer, :], axis=0)

    sel_layer = (np.abs(y_filtered) >= densf - delta_dens_flux) & (np.abs(y_filtered) < densf + delta_dens_flux)
    Flux_densf = np.mean(Flux_interp[sel_layer,:],axis=0)
    Flux_extended_densf = np.mean(Flux_extended_interp[sel_layer,:],axis=0)


    ########################################################################################################################
    # Here I calculate the carbon consumption rate due to (i) oxygen consumption and (ii) PARR
    ########################################################################################################################

    ############### I load Coriolis data with the oxygen information
    filename='6903095_Sprof_all.nc'
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
        date_time_obj = datetime.datetime.strptime(Date_Time_PARR[i], '%Y-%m-%dT%H:%M:%S')
        Date_Num_PARR[i] = calendar.timegm(date_time_obj.timetuple())
        #datetime.utcfromtimestamp(Date_Num[i])

    list_dates_PARR=np.unique(Date_Num_PARR)

    #############################################
    ############### Loop on the different isopycnal values chosen for the study of the oxygen profile
    Date_Num_limit=np.array([Date_Num.min(),Date_Num.min()+ndays]) #print(date_reference + datetime.timedelta(days=Date_Num.min()+127))

    # Here, for each profile included between two dates (Date_Num_limit), I compute the oxygen concentration between
    # dens0 and densf
    reference_isopycnal=(dens0+densf)*0.5
    reference_isopycnal_down=dens0
    reference_isopycnal_up=densf
    doxy_isopycnal=np.array([]);depth_isopycnal_tmp=np.array([]);depth_isopycnal_down_tmp=np.array([]);depth_isopycnal_up_tmp=np.array([])
    Date_Num_isopycnal=np.array([])
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
                depth_isopycnal_tmp2 = d[sel_layer]
                doxy_isopycnal = np.append(doxy_isopycnal, doxy_tmp)
                Date_Num_isopycnal = np.append(Date_Num_isopycnal, Date_Num[i])
                depth_isopycnal_tmp = np.append(depth_isopycnal_tmp, np.mean(depth_isopycnal_tmp2) )
                depth_isopycnal_down_tmp = np.append(depth_isopycnal_down_tmp, depth_isopycnal_tmp2[0] )
                depth_isopycnal_up_tmp = np.append(depth_isopycnal_up_tmp, depth_isopycnal_tmp2[-1] )
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
                    depth_isopycnal_down_tmp = np.append(depth_isopycnal_down_tmp, depth_isopycnal_tmp2[0])
                    depth_isopycnal_up_tmp = np.append(depth_isopycnal_up_tmp, depth_isopycnal_tmp2[-1])

    # If I have at least three points I interpolate them, so that the slope gives me the respiration rate
    # (in micromol/kg per day)
    depth_isopycnal = np.mean(depth_isopycnal_tmp)
    depth_isopycnal_down = np.mean(depth_isopycnal_down_tmp)
    depth_isopycnal_up = np.mean(depth_isopycnal_up_tmp)
    layer_thickness=np.mean(depth_isopycnal_up_tmp - depth_isopycnal_down_tmp)
    if Date_Num_isopycnal.size>2:
        (interpol,slpe_ci,_,signif,signif_label)=lin_fit(Date_Num_isopycnal,doxy_isopycnal)
        # fig = plt.figure(1, figsize=(12, 8))
        # plot1 = plt.scatter(Date_Num_isopycnal,doxy_isopycnal)
        slope_doxy = interpol.slope
        slope_ci_doxy = np.reshape(slpe_ci.copy(),(1,2))

        O2_resp_mgO_m3_d=-slope_doxy.copy()*reference_isopycnal*mol2gC/1000#*ndays*layer_thickness
        O2_resp_mgO_m3_d_ci=-slope_ci_doxy.copy()*( np.tile(reference_isopycnal,(2,1)).T) *mol2gC/1000#*ndays*layer_thickness
        O2_resp_mgC_m3_d = O2_resp_mgO_m3_d*Oxy2C
        #Since I have an uncertainty both on O2_resp_mgO_m3_d and on Oxy2C, I propagate the error to obtain the uncertainty on O2_resp_mgC_m3_d
        O2_resp_mgO_m3_d_ci=abs(np.diff(O2_resp_mgO_m3_d_ci)[0][0]/2)
        O2_resp_mgC_m3_d_err = np.sqrt( Oxy2C**2*O2_resp_mgO_m3_d_ci**2 + Oxy2C_std**2*O2_resp_mgO_m3_d**2 )
        O2_resp_mgC_m3_d_ci=np.ones((1,2));O2_resp_mgC_m3_d_ci[0,0]=O2_resp_mgC_m3_d-O2_resp_mgC_m3_d_err;O2_resp_mgC_m3_d_ci[0,1]=O2_resp_mgC_m3_d+O2_resp_mgC_m3_d_err

    else:
        O2_resp_mgC_m3_d=np.nan
        O2_resp_mgC_m3_d_ci=np.nan*np.ones((1,2))

    #############################################
    ############### Loop on different respiration types used to estimate PARR
    #List of the different Respiration types present in data
    list_Respi_types = [match for match in data.columns if "Respi" in match]
    nRespi= len(list_Respi_types)  #number of respiration types

    POC_resp_mgC_m3_d_list=np.zeros((nRespi,))
    POC_resp_mgC_m3_d_std_list = np.zeros((nRespi,))

    iRespi=0
    for iRespi in range(0,nRespi):
        PARR_nmol_l_h = np.array(data[list_Respi_types[iRespi]][sel_filename])
        # I convert the PARR measured in micromol/kg/day
        PARR_micromol_kg_day = PARR_nmol_l_h.copy() / 1000 * 24 / (dens_PARR / 1000)

        #############################################
        ############### Loop on the different PARR profiles

        #Here, for each profile included between two dates (Date_Num_limit), I compute PARR concentration between dens0
        # and densf
        PARR_isopycnal=np.array([]);#depth_PARR_tmp=np.array([]);dens_PARR_tmp=np.array([])
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
                    PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)
                else:  # If no values are found, then it could be that (if delta_rho is very small) the range reference_isopycnal_down–reference_isopycnal_up falls totally between two isopycnal layers: In that case, I extrapolate the PARR at that depth
                    for iy in range(0, y_PARR.size - 1):
                        if y_PARR[iy] <= reference_isopycnal < y_PARR[iy + 1]:
                            dist = (reference_isopycnal - y_PARR[iy]) / (y_PARR[iy + 1] - y_PARR[iy])
                            PARR_tmp = z_PARR[iy] + (z_PARR[iy + 1] - z_PARR[iy]) * dist
                            d_tmp = d_PARR[iy] + (d_PARR[iy + 1] - d_PARR[iy]) * dist
                            PARR_isopycnal = np.append(PARR_isopycnal, PARR_tmp)

        PARR_isopycnal_std = np.std(PARR_isopycnal)
        PARR_isopycnal = np.mean(PARR_isopycnal)

        # I convert the PARR and the oxygen respiration rates (in micromolO2/kg/d) to the total amount of carbon consumption
        # between depth0 and depthf, and between day0 and dayf (in mgC/m3/d)
        # *Oxy2C -> to micromolC/kg/d
        # *mol2gC -> to microgC/kg/d
        # /1000 -> to mgC/kg/d
        # *density -> to mgC/m3/d
        # *layer_thickness*ndays -> to mgC/m2

        # I calculate the PARR
        POC_resp_mgC_m3_d_list[iRespi] = PARR_isopycnal.copy()*reference_isopycnal*Oxy2C*mol2gC/1000 # *ndays * layer_thickness
        # I calculate the error on the PARR (POC_resp) by taking into account the error on Oxy2C as well
        POC_resp_mgC_m3_d_std_list[iRespi] = reference_isopycnal*mol2gC/1000*np.sqrt(PARR_isopycnal.copy()**2*Oxy2C_std**2+PARR_isopycnal_std.copy()**2*Oxy2C**2) # *ndays * layer_thickness

    # I convert bbp to mgC_m3_d: first, I calculate the error otherwise the formula is wrong
    # I calculate the error on the bbp_PARR by taking into account the error on Oxy2C as well
    bbp_PARR_dens0_densf_std = reference_isopycnal * mol2gC / 1000*np.sqrt(bbp_PARR_dens0_densf_std.copy()**2*Oxy2C**2 +bbp_PARR_dens0_densf.copy()**2*Oxy2C_std**2 ) # *ndays * layer_thickness
    bbp_PARR_Koestner_dens0_densf_std = reference_isopycnal * mol2gC / 1000*np.sqrt(bbp_PARR_Koestner_dens0_densf_std.copy()**2*Oxy2C**2 +bbp_PARR_Koestner_dens0_densf.copy()**2*Oxy2C_std**2 )   # *ndays * layer_thickness
    # I convert bbp to mgC_m3_d
    bbp_PARR_dens0_densf = bbp_PARR_dens0_densf.copy() * reference_isopycnal * Oxy2C * mol2gC / 1000  # *ndays * layer_thickness
    bbp_PARR_Koestner_dens0_densf = bbp_PARR_Koestner_dens0_densf.copy() * reference_isopycnal * Oxy2C * mol2gC / 1000  # *ndays * layer_thickness

    ########################################################################################################################
    # Here I calculate the carbon budget for depth0—depthf layer
    ########################################################################################################################
    # Date_Num_Flux = x_filtered
    # depth_POC_resp = list_depth_PARR
    # depth_02_resp = depth_isopycnal

    ############### I calculate the integrated POC (MiP+MaP+bbp), between depth0 and depthf, for day0 and dayf. I transform it to mgC/m2

    # t=np.r_[0:ndays]
    # (Integrated_POC_day0_mgC_m3_d,Integrated_POC_dayf_mgC_m3_d,Integrated_POC_day0_mgC_m3_d_std,Integrated_POC_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_mgC_m3)
    # (Integrated_POC_noBBP_day0_mgC_m3_d,Integrated_POC_noBBP_dayf_mgC_m3_d,Integrated_POC_noBBP_day0_mgC_m3_d_std,Integrated_POC_noBBP_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_noBBP_mgC_m3)
    # (Integrated_POC_Koestner_day0_mgC_m3_d,Integrated_POC_Koestner_dayf_mgC_m3_d,Integrated_POC_Koestner_day0_mgC_m3_d_std,Integrated_POC_Koestner_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_Koestner_mgC_m3)

    # I calculate the integrated POC between day and dayf, based on using the average value between day0 and day0+3days and between dayf-3days and dayf, respectively
    nd = 3
    Integrated_POC_day0_mgC_m3_d = np.nanmean(Integrated_POC_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_dayf_mgC_m3_d = np.nanmean(Integrated_POC_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_day0_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_noBBP_dayf_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_noBBP_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_day0_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_mgC_m3[0:nd]) / ndays# * layer_thickness
    Integrated_POC_Koestner_dayf_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_mgC_m3[-nd:]) / ndays# * layer_thickness
    Integrated_POC_Koestner_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Delta_Integrated_POC = Integrated_POC_dayf_mgC_m3_d - Integrated_POC_day0_mgC_m3_d
    Delta_Integrated_POC_std = np.sqrt( Integrated_POC_dayf_mgC_m3_d_std**2 + Integrated_POC_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_noBBP = Integrated_POC_noBBP_dayf_mgC_m3_d - Integrated_POC_noBBP_day0_mgC_m3_d
    Delta_Integrated_POC_noBBP_std = np.sqrt( Integrated_POC_noBBP_dayf_mgC_m3_d_std**2 + Integrated_POC_noBBP_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_Koestner = Integrated_POC_Koestner_dayf_mgC_m3_d - Integrated_POC_Koestner_day0_mgC_m3_d
    Delta_Integrated_POC_Koestner_std = np.sqrt( Integrated_POC_Koestner_dayf_mgC_m3_d_std**2 + Integrated_POC_Koestner_day0_mgC_m3_d_std**2 )

    # (Integrated_POC_extended_day0_mgC_m3_d,Integrated_POC_extended_dayf_mgC_m3_d,Integrated_POC_extended_day0_mgC_m3_d_std,Integrated_POC_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_extended_mgC_m3)
    # (Integrated_POC_noBBP_extended_day0_mgC_m3_d,Integrated_POC_noBBP_extended_dayf_mgC_m3_d,Integrated_POC_noBBP_extended_day0_mgC_m3_d_std,Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_noBBP_extended_mgC_m3)
    # (Integrated_POC_Koestner_extended_day0_mgC_m3_d,Integrated_POC_Koestner_extended_dayf_mgC_m3_d,Integrated_POC_Koestner_extended_day0_mgC_m3_d_std,Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std) = POC_day0dayf(t,Integrated_POC_Koestner_extended_mgC_m3)

    # Old way of calculate the integrated POC between day and dayf, based on using the exact value on day0 and dayf, respectively
    Integrated_POC_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_noBBP_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_day0_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_dayf_mgC_m3_d = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3[-nd:]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_day0_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3_std[0:nd]) / ndays # * layer_thickness
    Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std = np.nanmean(Integrated_POC_Koestner_extended_mgC_m3_std[-nd:]) / ndays # * layer_thickness
    Delta_Integrated_POC_extended = Integrated_POC_extended_dayf_mgC_m3_d - Integrated_POC_extended_day0_mgC_m3_d
    Delta_Integrated_POC_extended_std = np.sqrt( Integrated_POC_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_extended_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_noBBP_extended = Integrated_POC_noBBP_extended_dayf_mgC_m3_d - Integrated_POC_noBBP_extended_day0_mgC_m3_d
    Delta_Integrated_POC_noBBP_extended_std = np.sqrt( Integrated_POC_noBBP_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_noBBP_extended_day0_mgC_m3_d_std**2 )
    Delta_Integrated_POC_Koestner_extended = Integrated_POC_Koestner_extended_dayf_mgC_m3_d - Integrated_POC_Koestner_extended_day0_mgC_m3_d
    Delta_Integrated_POC_Koestner_extended_std = np.sqrt( Integrated_POC_Koestner_extended_dayf_mgC_m3_d_std**2 + Integrated_POC_Koestner_extended_day0_mgC_m3_d_std**2 )

    ############### I calculate the amount of POC entering from depht0 and exiting from dayf between day0 and dayf (in mgC/m3/day)

    # I extract the index of Flux_dens0/Flux_densf which correspond to day0 (and dayf)
    Flux_dens0_mgC_m3_d = np.mean(Flux_dens0) / layer_thickness # * ndays
    Flux_dens0_mgC_m3_d_std = np.std(Flux_dens0) / layer_thickness # * ndays
    Flux_densf_mgC_m3_d = np.mean(Flux_densf) / layer_thickness # * ndays
    Flux_densf_mgC_m3_d_std = np.std(Flux_densf) / layer_thickness # * ndays

    Flux_extended_dens0_mgC_m3_d = np.mean(Flux_extended_dens0) / layer_thickness # * ndays
    Flux_extended_dens0_mgC_m3_d_std = np.std(Flux_extended_dens0) / layer_thickness # * ndays
    Flux_extended_densf_mgC_m3_d = np.mean(Flux_extended_densf) / layer_thickness # * ndays
    Flux_extended_densf_mgC_m3_d_std = np.std(Flux_extended_densf) / layer_thickness # * ndays

    Delta_flux = Flux_dens0_mgC_m3_d - Flux_densf_mgC_m3_d
    Delta_flux_extended = Flux_extended_dens0_mgC_m3_d - Flux_extended_densf_mgC_m3_d

    Delta_flux_std = np.sqrt( Flux_dens0_mgC_m3_d_std**2 + Flux_densf_mgC_m3_d_std**2 )
    Delta_flux_extended_std = np.sqrt( Flux_extended_dens0_mgC_m3_d_std**2 + Flux_extended_densf_mgC_m3_d_std**2 )

    Theoretical_Budget = Delta_flux - Delta_Integrated_POC
    Theoretical_Budget_extended = Delta_flux_extended - Delta_Integrated_POC_extended
    Theoretical_Budget_noBBP = Delta_flux - Delta_Integrated_POC_noBBP
    Theoretical_Budget_noBBP_extended = Delta_flux_extended - Delta_Integrated_POC_noBBP_extended
    Theoretical_Budget_Koestner = Delta_flux - Delta_Integrated_POC_Koestner
    Theoretical_Budget_Koestner_extended = Delta_flux_extended - Delta_Integrated_POC_Koestner_extended

    Theoretical_Budget_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_std**2 )
    Theoretical_Budget_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_extended_std**2 )
    Theoretical_Budget_noBBP_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_noBBP_std**2 )
    Theoretical_Budget_noBBP_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_noBBP_extended_std**2 )
    Theoretical_Budget_Koestner_std = np.sqrt( Delta_flux_std**2 + Delta_Integrated_POC_Koestner_std**2 )
    Theoretical_Budget_Koestner_extended_std = np.sqrt( Delta_flux_extended_std**2 + Delta_Integrated_POC_Koestner_extended_std**2 )


    ############### I return the data
    return Theoretical_Budget,Theoretical_Budget_std, \
           Theoretical_Budget_extended,Theoretical_Budget_extended_std, \
           Theoretical_Budget_noBBP,Theoretical_Budget_noBBP_std, \
           Theoretical_Budget_noBBP_extended,Theoretical_Budget_noBBP_extended_std, \
           Theoretical_Budget_Koestner,Theoretical_Budget_Koestner_std, \
           Theoretical_Budget_Koestner_extended,Theoretical_Budget_Koestner_extended_std, \
           POC_resp_mgC_m3_d_list,POC_resp_mgC_m3_d_std_list,bbp_PARR_dens0_densf,bbp_PARR_dens0_densf_std,\
           bbp_PARR_Koestner_dens0_densf,bbp_PARR_Koestner_dens0_densf_std,O2_resp_mgC_m3_d,O2_resp_mgC_m3_d_ci,list_Respi_types,n_profiles, \
           Delta_Integrated_POC, Delta_Integrated_POC_std, Delta_Integrated_POC_noBBP, Delta_Integrated_POC_noBBP_std, Delta_Integrated_POC_Koestner, Delta_Integrated_POC_Koestner_std, \
           Delta_flux, Delta_flux_std,Flux_dens0_mgC_m3_d,Flux_densf_mgC_m3_d, \
           depth_isopycnal,depth_isopycnal_down,depth_isopycnal_up,layer_thickness,MiP_POC_dens0_densf,MiP_POC_extended_dens0_densf,MaP_POC_dens0_densf,bbp_POC_Koestner_dens0_densf
# endregion
#######################################################################
# Parameters for the carbon budget calculation
#######################################################################
day0=datetime.datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime.datetime(2021,7,31)        # starting date for the carbon budget calculation
ndays=(dayf-day0).days          # number of days
dens00=1026.3                   # starting isopycnal
dens_thickness=0.05             # thickness of the layer considered (in kg/m3)
delta_dens=0.025                 # every time I do a loop, how much I do increase depth0
densff=1027.5                   # final isopycnal investigated

dens0_list=np.r_[dens00:densff-dens_thickness+0.01:delta_dens]

#######################################################################
# I loop on the different depths
#######################################################################
Theoretical_Budget_list = np.array([])
Theoretical_Budget_extended_list = np.array([])
Theoretical_Budget_std_list = np.array([])
Theoretical_Budget_extended_std_list = np.array([])
Theoretical_Budget_noBBP_list = np.array([])
Theoretical_Budget_noBBP_extended_list = np.array([])
Theoretical_Budget_noBBP_std_list = np.array([])
Theoretical_Budget_noBBP_extended_std_list = np.array([])
Theoretical_Budget_Koestner_list = np.array([])
Theoretical_Budget_Koestner_extended_list = np.array([])
Theoretical_Budget_Koestner_std_list = np.array([])
Theoretical_Budget_Koestner_extended_std_list = np.array([])
POC_resp_mgC_m3_d_list = np.array([])
POC_resp_mgC_m3_d_std_list = np.array([])
bbpPARR_mgC_m3_d_list = np.array([])
bbpPARR_mgC_m3_d_std_list = np.array([])
bbpPARR_Koestner_mgC_m3_d_list = np.array([])
bbpPARR_Koestner_mgC_m3_d_std_list = np.array([])
O2_resp_mgC_m3_d_list = np.array([])
O2_resp_mgC_m3_d_ci_list = np.array([])
depth_isopycnal_list = np.array([])
depth_isopycnal_down_list = np.array([])
depth_isopycnal_up_list = np.array([])
layer_thickness_list = np.array([])
dens0=dens0_list[0]
for dens0 in dens0_list:
    densf = dens0 + dens_thickness
    (Theoretical_Budget,Theoretical_Budget_std,Theoretical_Budget_extended,Theoretical_Budget_extended_std,
       Theoretical_Budget_noBBP,Theoretical_Budget_noBBP_std,Theoretical_Budget_noBBP_extended,Theoretical_Budget_noBBP_extended_std,
       Theoretical_Budget_Koestner,Theoretical_Budget_Koestner_std,Theoretical_Budget_Koestner_extended,Theoretical_Budget_Koestner_extended_std,
       POC_resp_mgC_m3_d,POC_resp_mgC_m3_d_std,bbpPARR_mgC_m3_d,bbpPARR_mgC_m3_d_std,bbpPARR_Koestner_mgC_m3_d,bbpPARR_Koestner_mgC_m3_d_std,O2_resp_mgC_m3_d,O2_resp_mgC_m3_d_ci,RespirationTypes,n_profiles,
       _,_,_,_,_,_,_,_,_,_,depth_isopycnal,depth_isopycnal_down,depth_isopycnal_up,layer_thickness,_,_,_,_) = carbon_budget_calculation(dens0, densf, day0, dayf)

    Theoretical_Budget_list=np.append(Theoretical_Budget_list,Theoretical_Budget)
    Theoretical_Budget_extended_list=np.append(Theoretical_Budget_extended_list,Theoretical_Budget_extended)
    Theoretical_Budget_std_list=np.append(Theoretical_Budget_std_list,Theoretical_Budget_std)
    Theoretical_Budget_extended_std_list=np.append(Theoretical_Budget_extended_std_list,Theoretical_Budget_extended_std)
    Theoretical_Budget_noBBP_list=np.append(Theoretical_Budget_noBBP_list,Theoretical_Budget_noBBP)
    Theoretical_Budget_noBBP_extended_list=np.append(Theoretical_Budget_noBBP_extended_list,Theoretical_Budget_noBBP_extended)
    Theoretical_Budget_noBBP_std_list=np.append(Theoretical_Budget_noBBP_std_list,Theoretical_Budget_noBBP_std)
    Theoretical_Budget_noBBP_extended_std_list=np.append(Theoretical_Budget_noBBP_extended_std_list,Theoretical_Budget_noBBP_extended_std)
    Theoretical_Budget_Koestner_list=np.append(Theoretical_Budget_Koestner_list,Theoretical_Budget_Koestner)
    Theoretical_Budget_Koestner_extended_list=np.append(Theoretical_Budget_Koestner_extended_list,Theoretical_Budget_Koestner_extended)
    Theoretical_Budget_Koestner_std_list=np.append(Theoretical_Budget_Koestner_std_list,Theoretical_Budget_Koestner_std)
    Theoretical_Budget_Koestner_extended_std_list=np.append(Theoretical_Budget_Koestner_extended_std_list,Theoretical_Budget_Koestner_extended_std)
    POC_resp_mgC_m3_d_list=np.append(POC_resp_mgC_m3_d_list,POC_resp_mgC_m3_d,axis=0)
    POC_resp_mgC_m3_d_std_list=np.append(POC_resp_mgC_m3_d_std_list,POC_resp_mgC_m3_d_std,axis=0)
    bbpPARR_mgC_m3_d_list=np.append(bbpPARR_mgC_m3_d_list,bbpPARR_mgC_m3_d)
    bbpPARR_mgC_m3_d_std_list=np.append(bbpPARR_mgC_m3_d_std_list,bbpPARR_mgC_m3_d_std)
    bbpPARR_Koestner_mgC_m3_d_list=np.append(bbpPARR_Koestner_mgC_m3_d_list,bbpPARR_Koestner_mgC_m3_d)
    bbpPARR_Koestner_mgC_m3_d_std_list=np.append(bbpPARR_Koestner_mgC_m3_d_std_list,bbpPARR_Koestner_mgC_m3_d_std)
    O2_resp_mgC_m3_d_list=np.append(O2_resp_mgC_m3_d_list,O2_resp_mgC_m3_d)
    O2_resp_mgC_m3_d_ci_list=np.append(O2_resp_mgC_m3_d_ci_list,O2_resp_mgC_m3_d_ci.reshape((2,)),axis=0)
    depth_isopycnal_list=np.append(depth_isopycnal_list,depth_isopycnal)
    depth_isopycnal_down_list=np.append(depth_isopycnal_down_list,depth_isopycnal_down)
    depth_isopycnal_up_list=np.append(depth_isopycnal_up_list,depth_isopycnal_up)
    layer_thickness_list=np.append(layer_thickness_list,layer_thickness)

O2_resp_mgC_m3_d_ci_list=O2_resp_mgC_m3_d_ci_list.reshape(dens0_list.size,2)
POC_resp_mgC_m3_d_list=POC_resp_mgC_m3_d_list.reshape(dens0_list.size,len(RespirationTypes))
POC_resp_mgC_m3_d_std_list=POC_resp_mgC_m3_d_std_list.reshape(dens0_list.size,len(RespirationTypes))

dens_eddy_core_up = 1026.82
dens_eddy_core_down = 1027.2397618090454 #calculated at step 4 of Fig. 3a

#######################################################################
#region Plots
########################################################################################################################
######### Supplementary Fig. without extended size spectrum
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(2,3):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')
plt.plot(POC_resp_mgC_m3_d_list[:, 0] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 5] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_Koestner_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_list - Theoretical_Budget_Koestner_std_list, Theoretical_Budget_Koestner_list + Theoretical_Budget_Koestner_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,2)
# plt.ylabel('Dens (kg/m$^3$)', fontsize=fs)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
    yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an75/Q10_a_%0.2f_v07.pdf' % q10 ,dpi=200)
plt.close()


########################################################################################################################
######### Supplementary Fig. with extended size spectrum
########################################################################################################################
idx1,idx2=0,38
set_ylim_lower=depth_isopycnal_list[idx1]
set_ylim_upper=depth_isopycnal_list[idx2]
fig = plt.figure(2, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(O2_resp_mgC_m3_d_list,depth_isopycnal_list, 'k')
plt.scatter(O2_resp_mgC_m3_d_list,depth_isopycnal_list, c='black',s=5)
plt.fill_betweenx(depth_isopycnal_list, O2_resp_mgC_m3_d_ci_list[:, 1], O2_resp_mgC_m3_d_ci_list[:, 0], facecolor='b',color='gray', alpha=0.5, label='O$_2$ cons. rate')
for iResp in range(9,10):
    plt.plot(POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list, depth_isopycnal_list, c='b')

plt.fill_betweenx(depth_isopycnal_list, POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list-np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2),
                  POC_resp_mgC_m3_d_list[:,iResp] + bbpPARR_Koestner_mgC_m3_d_list+np.sqrt(POC_resp_mgC_m3_d_std_list[:,iResp]**2+bbpPARR_Koestner_mgC_m3_d_std_list**2), facecolor='b',
                  color='b', alpha=0.5, label='PARR\n($k_{rem}$=0.013d$^{-1}$;\nBelcher et al.)')

plt.plot(POC_resp_mgC_m3_d_list[:, 7] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='m',linestyle='dashed',label='PARR\n(Kalvelage\n/Iversen)')
plt.plot(POC_resp_mgC_m3_d_list[:, 12] + bbpPARR_Koestner_mgC_m3_d_list*0.1/0.013, depth_isopycnal_list, c='g',linestyle='dashed',label='PARR\n($k_{rem}$=0.1d$^{-1}$)')
plt.plot(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red')
plt.scatter(Theoretical_Budget_Koestner_extended_list, depth_isopycnal_list, c='red', s=5)
plt.fill_betweenx(depth_isopycnal_list, Theoretical_Budget_Koestner_extended_list - Theoretical_Budget_Koestner_extended_std_list, Theoretical_Budget_Koestner_extended_list + Theoretical_Budget_Koestner_extended_std_list,
                  facecolor='r', color='r', alpha=0.5, label='Bulk POC\nremov. rate')

plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_isopycnal_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.xlim(-0.05,2)
plt.xlabel('Carbon Consumption Rate (mgC/m$^3$/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
#I set yticks
nyticks=6
yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
yticks_down=np.linspace(depth_isopycnal_down_list[idx1], depth_isopycnal_down_list[idx2],nyticks)
yticks_up=np.linspace(depth_isopycnal_up_list[idx1], depth_isopycnal_up_list[idx2],nyticks)
yticklabels=[]
for i in range(0,nyticks):
        yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_isopycnal_list,dens0_list) ))
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(0.02, 1.075, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an75/Q10_b_%0.2f_v07.pdf' % q10 ,dpi=200)
plt.close()
#endregion
#endregion
