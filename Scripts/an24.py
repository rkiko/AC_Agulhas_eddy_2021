import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import cartopy
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename_mouw='%s/GIT/AC_Agulhas_eddy_2021/Data/an24/GO_flux.tab' % home

data=pd.read_csv(filename_mouw, sep='\t', header=87)
lon=data.Longitude
lat=data.Latitude
POC=data['POC flux [mg/m**2/day]']
Depth_trap=data['Depth water [m] (Sediment trap deployment depth)']
StartDate_trap=data['Date/Time (Deployed)']


########################################################################################################################
#Parameters to identify the area for which we want to extract the POC flux data
########################################################################################################################
lon0,lon1 = -10, 20
lat0, lat1 = -40,-20
Depth_max_trap = 60000
########################################################################################################################
#I extract the data according to the domain selection
########################################################################################################################
sel_in_domain=(lon>=lon0)&(lon<=lon1)&(lat>=lat0)&(lat<=lat1)&(Depth_trap<=Depth_max_trap)
print('Found %d traps' % sum(sel_in_domain))

POC_in_domain=POC[sel_in_domain]
lon_in_domain=lon[sel_in_domain]
lat_in_domain=lat[sel_in_domain]
Depth_trap_in_domain=Depth_trap[sel_in_domain]
StartDate_trap_in_domain = StartDate_trap[sel_in_domain]

print('Mean Flux is %0.1f mgC/m2/d' % np.mean(POC_in_domain))
print('Mean sediment trap depth is %0.0f' % np.mean(Depth_trap_in_domain))

width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(12, 8))
ax = fig.add_axes([0.12, 0.2, width, height])
plot = plt.scatter(POC_in_domain, Depth_trap_in_domain, c='blue')
plt.xlabel('POC Flux (mgC/m$^2$/d)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.ylim(0,4000),plt.xlim(0,np.max(POC_in_domain)*1.1)
plt.gca().invert_yaxis()
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/an24/POC_vs_Depth_an24.pdf', dpi=200)
plt.close()

delta_plot=1
set_xlim_lower, set_xlim_upper = lon0-delta_plot,lon1+delta_plot
set_ylim_lower, set_ylim_upper = lat0-delta_plot,lat1+delta_plot
fig = plt.figure(1, figsize=(12, 8))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
ax.set_extent([set_xlim_lower, set_xlim_upper, set_ylim_lower, set_ylim_upper])
ax.plot(lon_in_domain, lat_in_domain, 'k.', markersize=12)
ax.plot([5,18,18,5,5],[-36,-36,-32,-32,-36],'k',label='BGC Argo region')
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
gl.xlabel_style = {'fontsize': 15}
gl.ylabel_style = {'fontsize': 15}
ax.set_xlabel('Longitude', fontsize=15)
ax.set_ylabel('Latitude', fontsize=15)
gl.right_labels = False
gl.top_labels = False
ax.legend(fontsize=15)
plt.title('N traps: %d; Mean Flux: %0.1f mgC/m^2/d; Mean depth: %dm' % (sum(sel_in_domain),np.mean(POC_in_domain),np.mean(Depth_trap_in_domain)), fontsize=15)
plt.savefig('../Plots/an24/POC_trap_locations_an24.pdf', dpi=200)
plt.close()
