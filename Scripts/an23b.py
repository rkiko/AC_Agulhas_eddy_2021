import pickle
import numpy as np

from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

############### I load the variables

a_file = open("%s/an23/data_an23.pkl" % storedir, "rb")
data_an23 = pickle.load(a_file)
Integrated_POC_mgC_m3=data_an23['Integrated_POC_mgC_m3']
list_dates_Integrated_POC=data_an23['list_dates_Integrated_POC']
Flux_depth0=data_an23['Flux_depth0']
Flux_depthf=data_an23['Flux_depthf']
Date_Num_Flux=data_an23['Date_Num_Flux']
POC_resp_mgC_m2=data_an23['POC_resp_mgC_m2']
depth_POC_resp=data_an23['depth_POC_resp']
dens_POC_resp=data_an23['dens_POC_resp']
O2_resp_mgC_m2=data_an23['O2_resp_mgC_m2']
depth_02_resp=data_an23['depth_02_resp']
dens_02_resp=data_an23['dens_02_resp']
day0=data_an23['day0']
dayf=data_an23['dayf']
day0_float=data_an23['day0_float']
dayf_float=data_an23['dayf_float']
ndays=data_an23['ndays']
depth0=data_an23['depth0']
depthf=data_an23['depthf']
layer_thickness=data_an23['layer_thickness']
delta_depth=data_an23['delta_depth']
Oxy2C=data_an23['Oxy2C']
mol2gC=data_an23['mol2gC']


############### I calculate the integrated POC (MiP+MaP+bbp), between depth0 and depthf, for day0 and dayf. I transform it to mgC/m2

# I extract the index of Integrated_POC_mgC_m3 which correspond to day0 (and dayf)
tmp=list_dates_Integrated_POC-day0_float
idx0=np.where(np.abs(tmp)==(np.abs(tmp)).min())[0][0]
tmp=list_dates_Integrated_POC-dayf_float
idxf=np.where(np.abs(tmp)==(np.abs(tmp)).min())[0][0]

Integrated_POC_day0_mgC_m2 = Integrated_POC_mgC_m3[idx0]*layer_thickness
Integrated_POC_dayf_mgC_m2 = Integrated_POC_mgC_m3[idxf]*layer_thickness

############### I calculate the amount of POC entering from depht0 and exiting from dayf between day0 and dayf (in mgC/m2)

# I extract the index of Flux_depth0/Flux_depthf which correspond to day0 (and dayf)
tmp=Date_Num_Flux-day0_float
idx0=np.where(np.abs(tmp)==(np.abs(tmp)).min())[0][0]
tmp=Date_Num_Flux-dayf_float
idxf=np.where(np.abs(tmp)==(np.abs(tmp)).min())[0][0]

tmp=Flux_depth0[idx0:idxf]
Flux_depth0_mgC_m2=np.mean(tmp)*ndays
tmp=Flux_depthf[idx0:idxf]
Flux_depthf_mgC_m2=np.mean(tmp)*ndays

############### I calculate the amount of POC respired (in mgC/m2)
np.mean(O2_resp_mgC_m2)