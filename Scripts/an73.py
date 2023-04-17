
########################################################################################################################
########################################################################################################################
########################################################################################################################


# region an73
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

filename_eddy1Data = '%s/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv' % home
data_eddy1 = pd.read_csv(filename_eddy1Data, sep=',', header=0)
Date_Num_Eddy1 = data_eddy1['Datenum']  # mean radius_Vmax  = 60.25 km
radius_Vmax1 = data_eddy1['Rad_max']  # mean radius_Vmax  = 60.25 km
i1 = np.where(Date_Num_Eddy1 == matlab_datenum(2021, 4, 13))[0][0]
i2 = np.where(Date_Num_Eddy1 == matlab_datenum(2021, 7, 31))[0][0]
radius_mean = np.mean(radius_Vmax1[i1:i2 + 1]) * 1000  # in meters
radius_std = np.std(radius_Vmax1[i1:i2 + 1]) * 1000  # in meters
surface_eddy_mean = np.pi * radius_mean ** 2  # in square meters
surface_eddy_std = 2 * np.pi * radius_mean * radius_std  # in square meters

# I define the characteristics of our eddy (from an72) so that I don't have to calculate them again
POC_subducted_extended_0413to0731_450m_tonsC_day=219.33819491782526
subd_speed_eddyCore_0413to0731=1.099577867952028
Integrated_POC_Koestner_extended_mgC_m3=POC_subducted_extended_0413to0731_450m_tonsC_day/( subd_speed_eddyCore_0413to0731 * surface_eddy_mean)*10**9

#I define the characteristics of subsurface cylconic eddies in the 4 EBUS from Pegliasco et al. 2015. Thus, I use
#arrays of 4 elements. The first is the Benguela US, then there is the Canary Us, the California US, and finally the
#Peru Us. The data refer to long lived cyclones measured over 11 years
N_cyclones=np.array([474,396,511,461]) #total number of cyclones
frac_subsurface=np.array([0.55,0.6,0,0.12]) # fraction of the total number of cyclones which is subsurface intensified
radius_EBUS=np.array([100,95,83,98])*1000 #mean eddy radius of the cyclones in each EBUS (in meters)
radius_EBUS_std=np.array([46,37,38,41])*1000 #std of the mean eddy radius of the cyclones in each EBUS (in meters)
lifetime_EBUS=np.array([6.1,5.2,6.3,6.9])*30 #mean lifetime (in days) of the cyclones in each EBUS
lifetime_EBUS_std=np.array([7,6.3,7.5,7.5])*30 #std of the mean lifetime (in days) of the cyclones in each EBUS

#I calculate the number of subsurface intensified cyclones in each EBUS per year and their mean surface
N_subsurf_cyclones_perYear=N_cyclones*frac_subsurface/11
surface_US_mean = np.pi * radius_EBUS ** 2  # in square meters
surface_US_std = 2 * np.pi * radius_EBUS * radius_EBUS_std  # in square meters
#I calculate the mean export of POC of an average EBUS cyclone (in tonsC/day) assuming the same sinking speed than our
#cyclone and the same mean POC content
POC_subducted_EBUS_tonsC_day_by1cyclone = Integrated_POC_Koestner_extended_mgC_m3 * subd_speed_eddyCore_0413to0731 * surface_US_mean/10**9
#I calculate the same quantity but integrated over an year
POC_subducted_EBUS_tonsC_year_by1cyclone = POC_subducted_EBUS_tonsC_day_by1cyclone * lifetime_EBUS
#I sum up over all the subsurface intensifie cyclones
POC_subducted_EBUS_tonsC_year_tot = POC_subducted_EBUS_tonsC_year_by1cyclone * N_subsurf_cyclones_perYear
POC_subducted_EBUS_tonsC_year_totTOT=np.sum(POC_subducted_EBUS_tonsC_year_tot)
#I print the output
print("***********************************************************************************")
print("Estimation of full eddy core subduction pump (FECSP) in the 4 EBUS due to cyclones:")
print("***********************************************************************************");i=-1
i=i+1;print("FECSP in the Benguela upwelling system: %0.3e PgC/year" % (POC_subducted_EBUS_tonsC_year_tot[i]/10**9) )
i=i+1;print("FECSP in the Canary upwelling system: %0.3e PgC/year" % (POC_subducted_EBUS_tonsC_year_tot[i]/10**9) )
i=i+1;print("FECSP in the California upwelling system: %0.3e PgC/year" % (POC_subducted_EBUS_tonsC_year_tot[i]/10**9) )
i=i+1;print("FECSP in the Peru upwelling system: %0.3e PgC/year" % (POC_subducted_EBUS_tonsC_year_tot[i]/10**9) )
print("***********************************************************************************")
print("This results in a total of: %0.3e PgC/year by FECSP in all EBUS by cyclones" % (POC_subducted_EBUS_tonsC_year_totTOT/10**9) )
print("***********************************************************************************")
print("For comparison, the eddy subduction pump estimate of Boy et al (2019 ) is 0.09-2.0 PgC/year")
# endregion




