import os
import numpy as np
from scipy.special import lambertw
import pickle
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

#######################################################################
# In order to calculate the critical depth, we need the following quantities:
# 1) K, coefficient of attenuation
# 2) alpha^B
# 3) Irradiance along the whole day
# 4) Biomass Loss rate
# 5) W, the Lambert function
#######################################################################

#######################################################################
# I load
# 1) K, coefficient of attenuation
#######################################################################
a_file = open("%s/an42/data_an42.pkl" % storedir, "rb")
data_an42 = pickle.load(a_file)
Kd_490_float=data_an42['Kd_490_float']
Kd_490_datenum=data_an42['Kd_490_datenum']
a_file.close()
#I interpolate where I have a nan from neighbors
sel= np.where(np.isnan(Kd_490_float))
a = Kd_490_float[sel[0]-1]
b = Kd_490_float[sel[0]+1]
Kd_490_float[sel]=np.nanmean(np.array([a,b]),axis=0)
#######################################################################
# I load
# 2) alpha^B
#######################################################################
a_file = open("%s/an40/data_an40.pkl" % storedir, "rb")
data_an40 = pickle.load(a_file)
mean_alpha_in_domain=data_an40['mean_alpha_in_domain']
std_alpha_in_domain=data_an40['std_alpha_in_domain']
a_file.close()

#######################################################################
# I load
# 3) Irradiance along the whole day
#######################################################################
a_file = open("%s/an44/data_an44.pkl" % storedir, "rb")
data_an44 = pickle.load(a_file)
ssi_per_hour_float=data_an44['ssi_per_hour_float']
ssi_datenum=data_an44['ssi_datenum']
ssi_matlab_datenum=data_an44['ssi_matlab_datenum']
hour_daylight=data_an44['hour_daylight']
ssi_per_day_float = ssi_per_hour_float*hour_daylight
a_file.close()

#######################################################################
# I set
# 4) Biomass Loss rate
# According to Zhai et al. (2008, doi: 10.1029/2008GL035666)
#######################################################################
LTB=1.75*24 #mgC/(mg Chla)/day
LTB_std=2.5

#######################################################################
#######################################################################
# I calculate the critical depth according to the analytical formulation by Kovac et al. (2021, doi: 10.1093/icesjms/fsab013)
#######################################################################
#######################################################################
x = -mean_alpha_in_domain*ssi_per_day_float/LTB
x[x>-1]=np.nan
critical_depth = 1/Kd_490_float * (  lambertw( x*np.exp(x) ) - x  )
critical_depth = np.real(critical_depth)

x = -(mean_alpha_in_domain+std_alpha_in_domain*0.5)*ssi_per_day_float/(LTB-LTB_std*0.5)
x[x>-1]=np.nan
critical_depth_2 = 1/Kd_490_float * (  lambertw( x*np.exp(x) ) - x  )
critical_depth_2 = np.real(critical_depth_2)

x = -(mean_alpha_in_domain-std_alpha_in_domain*0.5)*ssi_per_day_float/(LTB+LTB_std*0.5)
x[x>-1]=np.nan
critical_depth_1 = 1/Kd_490_float * (  lambertw( x*np.exp(x) ) - x  )
critical_depth_1 = np.real(critical_depth_1)


dictionary_data = {"critical_depth": critical_depth,"critical_depth_1": critical_depth_1,"critical_depth_2": critical_depth_2,
                   "critical_depth_datenum": ssi_datenum,"critical_depth_matlab_datenum" : ssi_matlab_datenum}
a_file = open("%s/an45/data_an45.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()


#These lines are just to check how the critical depth changes when x increases (in absolute value). Answer: it increases
x=-np.r_[0.03:2.1:0.03];x=x[0:61]
critical_depth = 1/Kd_490_float * (  lambertw( x*np.exp(x) ) - x  )
critical_depth = np.real(critical_depth)
import matplotlib.pyplot as plt
plt.plot(x,critical_depth)
