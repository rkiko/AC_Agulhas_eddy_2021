
########################################################################################################################
########################################################################################################################
########################################################################################################################


# region an72
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
# I define the function for the carbon budget calculation
#######################################################################
# region carbon_subduction_calculation(dens0,densf,day0,dayf):
def carbon_subduction_calculation(dens0,densf,day0,dayf):
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

    # I convert the dates to float values (in seconds from 1970 1 1)
    Date_Num_bbp_calendar = Date_Num_bbp.copy()
    for i in range(0,Date_Num_bbp_calendar.size):
        date_time_obj = datetime.datetime(Date_Vec_bbp[i,0],Date_Vec_bbp[i,1],Date_Vec_bbp[i,2],
                                          Date_Vec_bbp[i,3],Date_Vec_bbp[i,4],Date_Vec_bbp[i,5])
        Date_Num_bbp_calendar[i] = calendar.timegm(date_time_obj.timetuple())
        # datetime.utcfromtimestamp(Date_Num[i])

    ########################################################################################################################
    # Here I calculate the mean integrated POC (i.e., MiP+MaP+bbp) in the isopycnal bin i.e. between dens0 and densf, and
    # averaged between day0 and dayf. To do so, (i) I filter it with a savgol function, then (ii) I interpolate it over a
    # regular grid versus time and density. This step is necessary to have MiP+MaP+bbp at 600 m, because some profiles only
    # reach 400 m; (iii) I extract the mean MiP+MaP+bbp values between dens0 and densf and between day0 and dayf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    MiP_filtered=np.array([]);dens_MiP_filtered=np.array([]);Date_Num_MiP_filtered=np.array([])
    MiP_extended_filtered=np.array([]);dens_MiP_extended_filtered=np.array([]);Date_Num_MiP_extended_filtered=np.array([])
    MaP_filtered=np.array([]);dens_MaP_filtered=np.array([]);Date_Num_MaP_filtered=np.array([])
    bbp_filtered=np.array([]);dens_bbp_filtered=np.array([]);Date_Num_bbp_filtered=np.array([])
    bbp_Koestner_filtered=np.array([]);dens_bbp_Koestner_filtered=np.array([]);Date_Num_bbp_Koestner_filtered=np.array([])

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


    ##############################################
    # Step 3, I calculate the mean MiP+MaP+bbp (and std) between dens0 and densf between day0 and dayf
    sel_dens0_densf = (np.abs(y_filtered) >= dens0) & (np.abs(y_filtered) < densf)
    MiP_POC_dens0_densf=np.mean(np.mean(MiP_interp[sel_dens0_densf,:],0))
    MiP_POC_extended_dens0_densf=np.mean(np.mean(MiP_extended_interp[sel_dens0_densf,:],0))
    MaP_POC_dens0_densf=np.mean(np.mean(MaP_interp[sel_dens0_densf,:],0))
    bbp_POC_dens0_densf=np.mean(np.mean(bbp_interp[sel_dens0_densf,:],0))
    bbp_POC_Koestner_dens0_densf=np.mean(np.mean(bbp_Koestner_interp[sel_dens0_densf,:],0))

    MiP_POC_dens0_densf_std = np.std(np.mean(MiP_interp[sel_dens0_densf,:],0))
    MiP_POC_extended_dens0_densf_std = np.std(np.mean(MiP_extended_interp[sel_dens0_densf,:],0))
    MaP_POC_dens0_densf_std = np.std(np.mean(MaP_interp[sel_dens0_densf,:],0))
    bbp_POC_dens0_densf_std = np.std(np.mean(bbp_interp[sel_dens0_densf,:],0))
    bbp_POC_Koestner_dens0_densf_std = np.std(np.mean(bbp_Koestner_interp[sel_dens0_densf,:],0))

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
    # Here I extract the average flux value in the isopycnal bin i.e. between dens0 and densf, and averaged between day0 and
    # dayf  at dens0 and densf. To do so, (i) I filter it with a savgol function, then (ii) I interpolate it over a regular
    # grid in time and density. This step is necessary to have the flux at 600 m, because some profiles only reach 400 m;
    # (iii) I extract the flux values between dens0 and densf and between day0 and dayf
    ########################################################################################################################

    ##############################################
    # Step 1 and 2, filter and interpolation
    Flux_filtered=np.array([]);dens_Flux_filtered=np.array([]);depth_Flux_filtered=np.array([]);Date_Num_Flux_filtered=np.array([])
    Flux_extended_filtered=np.array([]);dens_Flux_extended_filtered=np.array([]);Date_Num_Flux_extended_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i]
        z=Flux[sel];x=Date_Num[sel];y=dens[sel];d=depth[sel]
        sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2];d2=d[sel2]
        if sum(sel2) > 0:
            z = savgol_filter(z, 5, 1)
            Flux_filtered = np.concatenate((Flux_filtered, z))
            Date_Num_Flux_filtered = np.concatenate((Date_Num_Flux_filtered, x2))
            dens_Flux_filtered = np.concatenate((dens_Flux_filtered, y2))
            depth_Flux_filtered = np.concatenate((depth_Flux_filtered, d2))
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
    Depth_interp = griddata((Date_Num_Flux_extended_filtered, dens_Flux_extended_filtered), depth_Flux_filtered,(x_filtered_g, y_filtered_g), method="nearest")


    ##############################################
    # Step 3, flux extraction at dens0 and densf
    sel_dens0_densf = (np.abs(y_filtered) >= dens0) & (np.abs(y_filtered) < densf)
    Flux_mean = np.mean(np.mean(Flux_interp[sel_dens0_densf,:],axis=0))
    Flux_extended_mean = np.mean(np.mean(Flux_extended_interp[sel_dens0_densf, :], axis=0))
    Flux_std = np.std(np.mean(Flux_interp[sel_dens0_densf,:],axis=0))
    Flux_extended_std = np.std(np.mean(Flux_extended_interp[sel_dens0_densf, :], axis=0))

    ########################################################################################################################
    # Here I calculate the sinking speed of the isopycnal bin (in meters per day) and its mean width (in meters, anyway it's useless)
    ########################################################################################################################
    i=0;depth_vs_time=np.array([]);width_vs_time=np.array([])
    for i in range(0,x_filtered.size):
        sel_dens0_densf = (np.abs(y_filtered_g[:,i]) >= dens0) & (np.abs(y_filtered_g[:,i]) < densf)
        tmp = Depth_interp[sel_dens0_densf,i]
        depth_vs_time = np.append(depth_vs_time, np.mean(tmp) )
        width_vs_time = np.append(width_vs_time, tmp.max() - tmp.min() )

    (result, slope_ci, _, _, _) = lin_fit(x_filtered/86400,depth_vs_time)

    sink_speed=result.slope
    sink_speed_std=sink_speed-slope_ci[0]
    depth_mean = np.mean(depth_vs_time)
    width_mean = np.mean(width_vs_time)

    ########################################################################################################################
    # Here I calculate the mean eddy surface
    ########################################################################################################################
    filename_eddy1Data = '%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv' % home
    data_eddy1 = pd.read_csv(filename_eddy1Data, sep=',', header=0)
    Date_Num_Eddy1 = data_eddy1['Datenum']  # mean radius_Vmax  = 60.25 km
    radius_Vmax1 = data_eddy1['Rad_max']  # mean radius_Vmax  = 60.25 km
    i1 = np.where(Date_Num_Eddy1 == matlab_datenum(2021, 4, 13))[0][0]
    i2 = np.where(Date_Num_Eddy1 == matlab_datenum(2021, 7, 31))[0][0]
    radius_mean = np.mean(radius_Vmax1[i1:i2 + 1])*1000 #in meters
    radius_std  = np.std(radius_Vmax1[i1:i2 + 1])*1000 #in meters
    surface_eddy_mean=np.pi*radius_mean**2 #in square meters
    surface_eddy_std= 2*np.pi*radius_mean*radius_std #in square meters

    ########################################################################################################################
    # Here I calculate the POC subducted due to eddy subduction
    ########################################################################################################################
    POC_subducted_tonsC_day = Integrated_POC_Koestner_mgC_m3 * sink_speed * surface_eddy_mean/10**9
    POC_subducted_tonsC_day_std = np.sqrt(Integrated_POC_Koestner_mgC_m3_std**2 * sink_speed**2 * surface_eddy_mean**2 + Integrated_POC_Koestner_mgC_m3**2 * sink_speed_std**2 * surface_eddy_mean**2+ Integrated_POC_Koestner_mgC_m3**2 * sink_speed**2 * surface_eddy_std**2)/10**9
    POC_subducted_extended_tonsC_day = Integrated_POC_Koestner_extended_mgC_m3 * sink_speed * surface_eddy_mean/10**9
    POC_subducted_extended_tonsC_day_std = np.sqrt(Integrated_POC_Koestner_extended_mgC_m3_std**2 * sink_speed**2 * surface_eddy_mean**2 + Integrated_POC_Koestner_extended_mgC_m3**2 * sink_speed_std**2 * surface_eddy_mean**2+ Integrated_POC_Koestner_extended_mgC_m3**2 * sink_speed**2 * surface_eddy_std**2)/10**9
    POC_subducted_npBBP_tonsC_day = Integrated_POC_noBBP_mgC_m3 * sink_speed * surface_eddy_mean/10**9
    POC_subducted_npBBP_tonsC_day_std = np.sqrt(Integrated_POC_noBBP_mgC_m3_std**2 * sink_speed**2 * surface_eddy_mean**2 + Integrated_POC_noBBP_mgC_m3**2 * sink_speed_std**2 * surface_eddy_mean**2+ Integrated_POC_noBBP_mgC_m3**2 * sink_speed**2 * surface_eddy_std**2)/10**9

    ########################################################################################################################
    # Here I calculate the POC subducted due to gravitational pump (i.e. flux)
    ########################################################################################################################
    POC_BGP_tonsC_day = Flux_mean * surface_eddy_mean/10**9
    POC_BGP_tonsC_day_std = np.sqrt(Flux_std**2 * surface_eddy_mean**2 + Flux_mean**2 * surface_eddy_std**2)/10**9
    POC_BGP_extended_tonsC_day = Flux_extended_mean * surface_eddy_mean/10**9
    POC_BGP_extended_tonsC_day_std = np.sqrt(Flux_extended_std**2 * surface_eddy_mean**2 + Flux_extended_mean**2 * surface_eddy_std**2)/10**9
    ############### I return the data
    return POC_subducted_tonsC_day, POC_subducted_tonsC_day_std, POC_subducted_extended_tonsC_day, POC_subducted_extended_tonsC_day_std, \
           POC_subducted_npBBP_tonsC_day, POC_subducted_npBBP_tonsC_day_std, POC_BGP_tonsC_day, POC_BGP_tonsC_day_std, \
           POC_BGP_extended_tonsC_day, POC_BGP_extended_tonsC_day_std, sink_speed, sink_speed_std, depth_mean, width_mean
# endregion
#######################################################################
# Parameters for the carbon budget calculation
#######################################################################
day0=datetime.datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime.datetime(2021,7,31)        # starting date for the carbon budget calculation
ndays=(dayf-day0).days          # number of days
dens00=1026.35                   # starting isopycnal
dens_thickness=0.05             # thickness of the layer considered (in kg/m3)
delta_dens=0.025                 # every time I do a loop, how much I do increase depth0
densff=1027.35                   # final isopycnal investigated

dens0_list=np.r_[dens00:densff-dens_thickness+0.01:delta_dens]

#######################################################################
# I loop on the different depths
#######################################################################
POC_subducted_tonsC_day_list = np.array([])
POC_subducted_tonsC_day_std_list = np.array([])
POC_subducted_extended_tonsC_day_list = np.array([])
POC_subducted_extended_tonsC_day_std_list = np.array([])
POC_subducted_npBBP_tonsC_day_list = np.array([])
POC_subducted_npBBP_tonsC_day_std_list = np.array([])
POC_BGP_tonsC_day_list = np.array([])
POC_BGP_tonsC_day_std_list = np.array([])
POC_BGP_extended_tonsC_day_list = np.array([])
POC_BGP_extended_tonsC_day_std_list = np.array([])
sink_speed_list = np.array([])
sink_speed_std_list = np.array([])
depth_list = np.array([])
width_list = np.array([])
dens0=dens0_list[0]
ct=0
for dens0 in dens0_list:
    ct=ct+1
    print('Loop %d out of %d' %(ct,dens0_list.size))
    densf = dens0 + dens_thickness
    (POC_subducted_tonsC_day, POC_subducted_tonsC_day_std, POC_subducted_extended_tonsC_day, POC_subducted_extended_tonsC_day_std,
     POC_subducted_npBBP_tonsC_day, POC_subducted_npBBP_tonsC_day_std, POC_BGP_tonsC_day, POC_BGP_tonsC_day_std,
     POC_BGP_extended_tonsC_day, POC_BGP_extended_tonsC_day_std, sink_speed, sink_speed_std, depth_mean, width_mean) = carbon_subduction_calculation(dens0, densf, day0, dayf)

    POC_subducted_tonsC_day_list=np.append(POC_subducted_tonsC_day_list,POC_subducted_tonsC_day)
    POC_subducted_tonsC_day_std_list=np.append(POC_subducted_tonsC_day_std_list,POC_subducted_tonsC_day_std)
    POC_subducted_extended_tonsC_day_list=np.append(POC_subducted_extended_tonsC_day_list,POC_subducted_extended_tonsC_day)
    POC_subducted_extended_tonsC_day_std_list=np.append(POC_subducted_extended_tonsC_day_std_list,POC_subducted_extended_tonsC_day_std)
    POC_subducted_npBBP_tonsC_day_list=np.append(POC_subducted_npBBP_tonsC_day_list,POC_subducted_npBBP_tonsC_day)
    POC_subducted_npBBP_tonsC_day_std_list=np.append(POC_subducted_npBBP_tonsC_day_std_list,POC_subducted_npBBP_tonsC_day_std)
    POC_BGP_tonsC_day_list=np.append(POC_BGP_tonsC_day_list,POC_BGP_tonsC_day)
    POC_BGP_tonsC_day_std_list=np.append(POC_BGP_tonsC_day_std_list,POC_BGP_tonsC_day_std)
    POC_BGP_extended_tonsC_day_list=np.append(POC_BGP_extended_tonsC_day_list,POC_BGP_extended_tonsC_day)
    POC_BGP_extended_tonsC_day_std_list=np.append(POC_BGP_extended_tonsC_day_std_list,POC_BGP_extended_tonsC_day_std)
    sink_speed_list=np.append(sink_speed_list,sink_speed)
    sink_speed_std_list=np.append(sink_speed_std_list,sink_speed_std)
    depth_list=np.append(depth_list,depth_mean)
    width_list=np.append(width_list,width_mean)

dens_eddy_core_up = 1026.82
dens_eddy_core_down = 1027.2397618090454 #calculated at step 4 of Fig. 3a

#######################################################################
# I extract the flux at 100m to calculate the Martin curve which I add to the plot for comparison
#######################################################################

(POC_subducted_tonsC_day, POC_subducted_tonsC_day_std, POC_subducted_extended_tonsC_day, POC_subducted_extended_tonsC_day_std,
 POC_subducted_npBBP_tonsC_day, POC_subducted_npBBP_tonsC_day_std, POC_BGP_tonsC_day, POC_BGP_tonsC_day_std,
 POC_BGP_extended_tonsC_day, POC_BGP_extended_tonsC_day_std, sink_speed, sink_speed_std, depth_mean, width_mean) = carbon_subduction_calculation(dens0, densf, day0, dayf)

#######################################################################
#region Plots
########################################################################################################################
######### Fig.: with BBP from Koestner
########################################################################################################################
idx1,idx2=0,depth_list.size-1
set_ylim_lower=depth_list[idx1]
set_ylim_upper=depth_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(POC_subducted_tonsC_day_list,depth_list, 'r')
plt.scatter(POC_subducted_tonsC_day_list,depth_list, c='red',s=5)
plt.fill_betweenx(depth_list, POC_subducted_tonsC_day_list-POC_subducted_tonsC_day_std_list, POC_subducted_tonsC_day_list+POC_subducted_tonsC_day_std_list, facecolor='b',color='red', alpha=0.5, label='FESP')
plt.plot(POC_BGP_tonsC_day_list,depth_list, 'b')
plt.scatter(POC_BGP_tonsC_day_list,depth_list, c='b',s=5)
plt.fill_betweenx(depth_list, POC_BGP_tonsC_day_list-POC_BGP_tonsC_day_std_list, POC_BGP_tonsC_day_list+POC_BGP_tonsC_day_std_list, facecolor='b',color='blue', alpha=0.5, label='BGP')
plt.xlim(0,1000)
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.ylabel('Depth (m)', fontsize=fs)
plt.xlabel('POC export rate (tons C/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
# #I set yticks
# nyticks=6
# yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
# yticks_down=np.linspace((depth_list-width_list*0.5)[idx1], (depth_list-width_list*0.5)[idx2],nyticks)
# yticks_up=np.linspace((depth_list+width_list*0.5)[idx1], (depth_list+width_list*0.5)[idx2],nyticks)
# yticklabels=[]
# for i in range(0,nyticks):
#     yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_list,dens0_list) ))
# ax.set_yticks(yticks)
# ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
# ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an72/POC_sequestered_an72.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Fig.: Extended size spectrum with BBP from Koestner
########################################################################################################################
idx1,idx2=0,depth_list.size-1
set_ylim_lower=depth_list[idx1]
set_ylim_upper=depth_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(POC_subducted_extended_tonsC_day_list,depth_list, 'r')
plt.scatter(POC_subducted_extended_tonsC_day_list,depth_list, c='red',s=5)
plt.fill_betweenx(depth_list, POC_subducted_extended_tonsC_day_list-POC_subducted_extended_tonsC_day_std_list, POC_subducted_extended_tonsC_day_list+POC_subducted_extended_tonsC_day_std_list, facecolor='b',color='red', alpha=0.5, label='FESP')
plt.plot(POC_BGP_extended_tonsC_day_list,depth_list, 'b')
plt.scatter(POC_BGP_extended_tonsC_day_list,depth_list, c='b',s=5)
plt.fill_betweenx(depth_list, POC_BGP_extended_tonsC_day_list-POC_BGP_extended_tonsC_day_std_list, POC_BGP_extended_tonsC_day_list+POC_BGP_extended_tonsC_day_std_list, facecolor='b',color='blue', alpha=0.5, label='BGP')
plt.xlim(0,1000)
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(depth_list[1], xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod',linestyles='dotted',linewidth=5,zorder=20)
plt.ylabel('Depth (m)', fontsize=fs)
plt.xlabel('POC export rate (tons C/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
# #I set yticks
# nyticks=6
# yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
# yticks_down=np.linspace((depth_list-width_list*0.5)[idx1], (depth_list-width_list*0.5)[idx2],nyticks)
# yticks_up=np.linspace((depth_list+width_list*0.5)[idx1], (depth_list+width_list*0.5)[idx2],nyticks)
# yticklabels=[]
# for i in range(0,nyticks):
#     yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_list,dens0_list) ))
# ax.set_yticks(yticks)
# ax.set_yticklabels(yticklabels,fontsize=6)
ax.text(-0.15, 1.005, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
# ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an72/POC_sequestered_extended_an72.pdf' ,dpi=200)
plt.close()


########################################################################################################################
######### Supplementary Fig. 04A: without BBP
########################################################################################################################
idx1,idx2=0,depth_list.size-1
set_ylim_lower=depth_list[idx1]
set_ylim_upper=depth_list[idx2]
fs=10
width, height = 0.72, 0.8
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.23, 0.15, width, height], ylim=(set_ylim_lower, set_ylim_upper))
plt.plot(POC_subducted_npBBP_tonsC_day_list,depth_list, 'r')
plt.scatter(POC_subducted_npBBP_tonsC_day_list,depth_list, c='red',s=5)
plt.fill_betweenx(depth_list, POC_subducted_npBBP_tonsC_day_list-POC_subducted_npBBP_tonsC_day_std_list, POC_subducted_npBBP_tonsC_day_list+POC_subducted_npBBP_tonsC_day_std_list, facecolor='b',color='red', alpha=0.5, label='Subducted')
plt.plot(POC_BGP_tonsC_day_list,depth_list, 'b')
plt.scatter(POC_BGP_tonsC_day_list,depth_list, c='b',s=5)
plt.fill_betweenx(depth_list, POC_BGP_tonsC_day_list-POC_BGP_tonsC_day_std_list, POC_BGP_tonsC_day_list+POC_BGP_tonsC_day_std_list, facecolor='b',color='blue', alpha=0.5, label='Gravitational pump')
plt.xlim(0,1000)
plt.hlines(200, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
plt.hlines(600, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color='darkgoldenrod')
# plt.ylabel('Dens (kg/m$^3$)', fontsize=fs)
plt.xlabel('POC export rate (tons C/d)', fontsize=fs)
plt.legend(fontsize=7)
plt.gca().invert_yaxis()
# #I set yticks
# nyticks=6
# yticks=np.linspace(set_ylim_lower, set_ylim_upper,nyticks)
# yticks_down=np.linspace((depth_list-width_list*0.5)[idx1], (depth_list-width_list*0.5)[idx2],nyticks)
# yticks_up=np.linspace((depth_list+width_list*0.5)[idx1], (depth_list+width_list*0.5)[idx2],nyticks)
# yticklabels=[]
# for i in range(0,nyticks):
#     yticklabels.append('[%d–%dm]\n%0.2f kg/m$^3$' % (yticks_down[i],yticks_up[i], np.interp(yticks[i],depth_list,dens0_list) ))
# ax.set_yticks(yticks)
# ax.set_yticklabels(yticklabels,fontsize=6)
# ax.text(-0.25, 1.075, 'a', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
# ax.text(1.075, 1.06, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an72/POC_sequestered_noBBP_an72.pdf' ,dpi=200)
plt.close()



#endregion
#######################################################################

#######################################################################
#region Write some values for latex doc
#######################################################################

from write_latex_data import write_latex_data
filename='%s/GIT/AC_Agulhas_eddy_2021/Data/data_latex_Agulhas.dat' % home
#These values are extracted from the bulk POC and PARR calculated in the eddy core considering only one unique layer
depth2 = np.r_[depth_list.min():depth_list.max():10]
sel_inEddyCore = (depth2>200)&(depth2<600)
tmp = np.interp(depth2,depth_list,POC_BGP_extended_tonsC_day_list)
argument = 'BGP_extended_0413to0731_tonsC_day'
BGP_extended_0413to0731_tonsC_day=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%d' % BGP_extended_0413to0731_tonsC_day)
tmp = np.interp(depth2,depth_list,POC_BGP_extended_tonsC_day_std_list)
argument = 'BGP_extended_std_0413to0731_tonsC_day'
BGP_extended_std_0413to0731_tonsC_day=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%d' % BGP_extended_std_0413to0731_tonsC_day)
tmp = np.interp(depth2,depth_list,POC_subducted_extended_tonsC_day_list)
argument = 'POC_subducted_extended_0413to0731_tonsC_day'
POC_subducted_extended_0413to0731_tonsC_day=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%d' % POC_subducted_extended_0413to0731_tonsC_day)
tmp = np.interp(depth2,depth_list,POC_subducted_extended_tonsC_day_std_list)
argument = 'POC_subducted_extended_std_0413to0731_tonsC_day'
POC_subducted_extended_std_0413to0731_tonsC_day=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%d' % POC_subducted_extended_std_0413to0731_tonsC_day)
argument = 'Tot_POC_export_extended_0413to0731_tonsC_day'
Tot_POC_export_extended_0413to0731_tonsC_day= BGP_extended_0413to0731_tonsC_day + POC_subducted_extended_0413to0731_tonsC_day
write_latex_data(filename,argument,'%d' % Tot_POC_export_extended_0413to0731_tonsC_day)
argument = 'Tot_POC_export_extended_std_0413to0731_tonsC_day'
arg_value= np.sqrt(BGP_extended_std_0413to0731_tonsC_day**2 + POC_subducted_extended_std_0413to0731_tonsC_day**2)
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'Ratio_POC_subduction2BGP_extended_0413to0731_tonsC_day'
arg_value= POC_subducted_extended_0413to0731_tonsC_day/BGP_extended_0413to0731_tonsC_day*100
write_latex_data(filename,argument,'%d' % arg_value)
argument = 'nPeople_equivalentCO2emitted'
arg_value= Tot_POC_export_extended_0413to0731_tonsC_day/(6.9/365)
write_latex_data(filename,argument,'%d' % arg_value)

tmp = np.interp(depth2,depth_list,sink_speed_list)
argument = 'subd_speed_eddyCore_0413to0731'
arg_value=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%0.1f' % arg_value)
tmp = np.interp(depth2,depth_list,sink_speed_std_list)
argument = 'subd_speed_eddyCore_std_0413to0731'
arg_value=np.mean(tmp[sel_inEddyCore])
write_latex_data(filename,argument,'%0.1f' % arg_value)

depth2write = np.array([200,450])
depth_tmp = 200
for depth_tmp in depth2write:
    argument = 'BGP_extended_0413to0731_%dm_tonsC_day' % (depth_tmp)
    BGP = np.interp(depth_tmp,depth_list,POC_BGP_extended_tonsC_day_list)
    write_latex_data(filename,argument,'%d' %  BGP)
    argument = 'BGP_extended_std_0413to0731_%dm_tonsC_day' % (depth_tmp)
    BGP_std = np.interp(depth_tmp,depth_list,POC_BGP_extended_tonsC_day_std_list)
    write_latex_data(filename,argument,'%d' % BGP_std )
    argument = 'POC_subducted_extended_0413to0731_%dm_tonsC_day' % (depth_tmp)
    FESP = np.interp(depth_tmp,depth_list,POC_subducted_extended_tonsC_day_list)
    write_latex_data(filename,argument,'%d' %  FESP)
    argument = 'POC_subducted_extended_std_0413to0731_%dm_tonsC_day' % (depth_tmp)
    FESP_std = np.interp(depth_tmp,depth_list,POC_subducted_extended_tonsC_day_std_list)
    write_latex_data(filename,argument,'%d' %  FESP_std)
    argument = 'Tot_POC_export_extended_0413to0731_%dm_tonsC_day' % (depth_tmp)
    write_latex_data(filename,argument,'%d' %  np.round(FESP+BGP))
    argument = 'Tot_POC_export_extended_std_0413to0731_%dm_tonsC_day' % (depth_tmp)
    write_latex_data(filename,argument,'%d' %  np.round(np.sqrt(FESP_std**2+BGP_std**2)))
    argument = 'Ratio_POC_subduction2BGP_extended_0413to0731_%dm_tonsC_day' % (depth_tmp)
    write_latex_data(filename, argument, '%d' % np.round(FESP/BGP*100))

# endregion
#######################################################################
# endregion



