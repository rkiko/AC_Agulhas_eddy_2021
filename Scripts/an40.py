import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import cartopy
import pickle
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename_alpha='%s/GIT/AC_Agulhas_eddy_2021/Data/an40/Bouman_2017.tab' % home

data=pd.read_csv(filename_alpha, sep='\t', header=57)
lon=data.Longitude
lat=data.Latitude
alpha=data['alpha [(mg C/mg Chl a/h)/(µE/m**2/s)]']
depth=data['Depth water [m]']
Date_=data['Date/Time']
chief_scientist=data['Chief scientist(s)']

#I plot the stations
# delta_plot=1
set_xlim_lower, set_xlim_upper = -20,30
set_ylim_lower, set_ylim_upper = -45,-25
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
# ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.plot(lon, lat, 'k.', markersize=12)
ax.plot([5,18,18,5,5],[-36,-36,-32,-32,-36],'k',label='BGC Argo region')
gl.xlabel_style = {'fontsize': 15}
gl.ylabel_style = {'fontsize': 15}
ax.set_xlabel('Longitude', fontsize=15)
ax.set_ylabel('Latitude', fontsize=15)
gl.right_labels = False
gl.top_labels = False
ax.legend(fontsize=15)
# plt.title('N traps: %d; Mean Flux: %0.1f mgC/m^2/d; Mean depth: %dm' % (sum(sel_in_domain),np.mean(POC_in_domain),np.mean(Depth_trap_in_domain)), fontsize=15)
plt.savefig('../Plots/an40/Alpha_station_locations_global_an40.pdf', dpi=200)
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
plt.savefig('../Plots/an40/Alpha_station_locations_zoom_an40.pdf', dpi=200)

plt.close()


########################################################################################################################
#Parameters to identify the area for which we want to extract Alpha^B
########################################################################################################################

lon0,lon1 = -10, 20
lat0, lat1 = -39,-24
########################################################################################################################
#I extract the data according to the domain selection
########################################################################################################################
sel_in_domain=(lon>=lon0)&(lon<=lon1)&(lat>=lat0)&(lat<=lat1)
print('Found %d stations' % sum(sel_in_domain))

alpha_in_domain=alpha[sel_in_domain]
lon_in_domain=lon[sel_in_domain]
lat_in_domain=lat[sel_in_domain]
Depth_in_domain=depth[sel_in_domain]
Date_in_domain = Date_[sel_in_domain]

mean_alpha_in_domain = np.mean(alpha_in_domain)
std_alpha_in_domain = np.std(alpha_in_domain)
print('Mean alpha is %0.3f mgC/(mg Chla)/h' % np.mean(alpha_in_domain))
print('Mean station depth is %0.0f' % np.mean(Depth_in_domain))

storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home

dictionary_data = {"mean_alpha_in_domain": mean_alpha_in_domain, "std_alpha_in_domain": std_alpha_in_domain}
a_file = open("%s/an40/data_an40.pkl" % storedir, "wb")
pickle.dump(dictionary_data, a_file)
a_file.close()
