import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from datetime import datetime
import calendar
from scipy.signal import savgol_filter
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
########################################################################################################################
# Flux from Mauw et al.
########################################################################################################################
########################################################################################################################

########################################################################################################################
#Parameters to identify the area for which we want to extract the POC flux data
########################################################################################################################
lon0,lon1 = -10, 20
lat0, lat1 = -40,-20
Depth_max_trap = 600
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


########################################################################################################################
########################################################################################################################
# Flux from our BGC Argo float
########################################################################################################################
########################################################################################################################

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
day0=datetime(2021,4,13)        # starting date for the carbon budget calculation
dayf=datetime(2021,10,18)       # final date for the carbon budget calculation
depthf=600                      # final depth
delta_depth=15                  # around of the depth which I consider when extracting the flux
day0_float=calendar.timegm(day0.timetuple())
dayf_float=calendar.timegm(dayf.timetuple())

########################################################################################################################
# I process the data
########################################################################################################################
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
Flux = np.array(data['Flux_mgC_m2'][sel_filename])
Flux_Clements=np.array(data['Flux_Clements_mgC_m2'][sel_filename])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_filename])
Flux_extended = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_filename])
Flux_extended_eta_b = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_filename])
Flux_Clements_extended=np.array(data['Flux_Clements_mgC_m2_from0.0254sizeclass'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I select the data only in the prescribed period
list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]

##############################################
# I filter the Flux
Flux_filtered_depthf=np.array([])
Flux_Clements_filtered_depthf=np.array([])
Flux_eta_b_filtered_depthf=np.array([])
Flux_extended_filtered_depthf=np.array([])
Flux_extended_eta_b_filtered_depthf=np.array([])
Flux_Clements_extended_filtered_depthf=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_filtered_depthf = np.append(Flux_filtered_depthf, np.mean(z[sel_depthf]) )
    z=Flux_Clements[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_Clements_filtered_depthf = np.append(Flux_Clements_filtered_depthf, np.mean(z[sel_depthf]) )
    z=Flux_eta_b[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_eta_b_filtered_depthf = np.append(Flux_eta_b_filtered_depthf, np.mean(z[sel_depthf]) )
    z=Flux_extended[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_extended_filtered_depthf = np.append(Flux_extended_filtered_depthf, np.mean(z[sel_depthf]) )
    z=Flux_Clements_extended[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_Clements_extended_filtered_depthf = np.append(Flux_Clements_extended_filtered_depthf, np.mean(z[sel_depthf]) )
    z=Flux_extended_eta_b[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_extended_eta_b_filtered_depthf = np.append(Flux_extended_eta_b_filtered_depthf, np.mean(z[sel_depthf]) )



########################################################################################################################
########################################################################################################################
# Flux from South Africa - South America transect
########################################################################################################################
########################################################################################################################

########################################################################################################################
#Parameters for the carbon budget calculation
########################################################################################################################
lon0=10                          # Longitudinal limits of
lonf=18                         # the flux data
depthf=600                      # Final depth
delta_depth=15                  # Around of the depth which I consider when extracting the flux

########################################################################################################################
# I process the data
########################################################################################################################
filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_diagnostics_data_231.tsv' % home
data=pd.read_csv(filename_ecopart, sep='\t', header=0)
Lon=np.array(data['Longitude'])

#I select only the profiles which were taken in the longitudal band I wish
sel_longitude = (Lon>lon0)&(Lon<lonf)

# I extract the data
Date_Time=np.array(data['Date_Time'][sel_longitude])
depth=np.array(data['Depth [m]'][sel_longitude])
Flux = np.array(data['Flux_mgC_m2'][sel_longitude])
Flux_eta_b=np.array(data['Flux_mgC_m2_from0.1200sizeclass_eta0.62_b66'][sel_longitude])
Flux_extended = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b132'][sel_longitude])
Flux_extended_eta_b = np.array(data['Flux_mgC_m2_from0.0254sizeclass_eta0.62_b66'][sel_longitude])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I select the data only in the prescribed period
list_dates=np.sort(np.unique(Date_Num))

##############################################
# I filter the Flux
Flux_filtered_depthf_231=np.array([])
Flux_eta_b_filtered_depthf_231=np.array([])
Flux_extended_filtered_depthf_231=np.array([])
Flux_extended_eta_b_filtered_depthf_231=np.array([])
i=0
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i]
    z=Flux[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_filtered_depthf_231 = np.append(Flux_filtered_depthf_231, np.mean(z[sel_depthf]) )
    z=Flux_eta_b[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_eta_b_filtered_depthf_231 = np.append(Flux_eta_b_filtered_depthf_231, np.mean(z[sel_depthf]) )
    z=Flux_extended[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_extended_filtered_depthf_231 = np.append(Flux_extended_filtered_depthf_231, np.mean(z[sel_depthf]) )
    z=Flux_extended_eta_b[sel];x=Date_Num[sel];y=depth[sel]
    sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if (sum(sel2) > 0)&(y.max()>depthf):
        z = savgol_filter(z, 5, 1)
        sel_depthf=(np.abs(y)>depthf-delta_depth)&(np.abs(y)<=depthf+delta_depth)
        Flux_extended_eta_b_filtered_depthf_231 = np.append(Flux_extended_eta_b_filtered_depthf_231, np.mean(z[sel_depthf]) )




########################################################################################################################
########################################################################################################################
# I plot flux from our data and from literature
########################################################################################################################
########################################################################################################################
fs=8
width, height = 0.78, 0.72

# First plot: flux calculated without considering smallest size classes and with old eta and b values
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.18, width, height])
plt.boxplot([Flux_filtered_depthf,Flux_Clements_filtered_depthf,Flux_filtered_depthf_231,POC_in_domain])
plt.ylim(0,70)
plt.ylabel('POC Flux (mgC/m$^2/d$)', fontsize=fs)
plt.title('Argo Flux , eta=0.62,b=132\nbetween %d-%02d-%02d and %d-%02d-%02d' % (day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day), fontsize=8)
plt.xticks([1,2,3,4],['BGC Argo\n6903095 ','BGC Argo\n6903095\n(Clements\net al.)','Transect\nMSM60','Sediment traps\n(Mouw et al.)'], fontsize=fs)
plt.savefig('../Plots/an32/01oldway_POCFlux_Argo_vs_literature_%d%02d%02dto%d%02d%02d_an32.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Second plot: flux calculated without considering smallest size classes and with new eta and b values
fig = plt.figure(2, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.18, width, height])
plt.boxplot([Flux_eta_b_filtered_depthf,Flux_Clements_filtered_depthf,Flux_eta_b_filtered_depthf_231,POC_in_domain])
plt.ylim(0,70)
plt.ylabel('POC Flux (mgC/m$^2/d$)', fontsize=fs)
plt.title('Argo Flux no small size classes, eta=0.62,b=66\n between %d-%02d-%02d and %d-%02d-%02d' % (day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day), fontsize=8)
plt.xticks([1,2,3,4],['BGC Argo\n6903095 ','BGC Argo\n6903095\n(Clements\net al.)','Transect\nMSM60','Sediment traps\n(Mouw et al.)'], fontsize=fs)
ax.text(-0.15, 1.125, 'a', transform=ax.transAxes, fontsize=12, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.savefig('../Plots/an32/02eta_b_POCFlux_Argo_vs_literature_%d%02d%02dto%d%02d%02d_an32.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Third plot: flux calculated considering smallest size classes but with old eta and b values
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.18, width, height])
plt.boxplot([Flux_extended_filtered_depthf,Flux_Clements_extended_filtered_depthf,Flux_extended_filtered_depthf_231,POC_in_domain])
plt.ylim(0,70)
plt.ylabel('POC Flux (mgC/m$^2/d$)', fontsize=fs)
plt.title('Argo Flux with small size classes, eta=0.62,b=132\n between %d-%02d-%02d and %d-%02d-%02d' % (day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day), fontsize=8)
plt.xticks([1,2,3,4],['BGC Argo\n6903095 ','BGC Argo\n6903095\n(Clements\net al.)','Transect\nMSM60','Sediment traps\n(Mouw et al.)'], fontsize=fs)
plt.savefig('../Plots/an32/03extended_POCFlux_Argo_vs_literature_%d%02d%02dto%d%02d%02d_an32.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()

# Fourth plot: flux calculated considering smallest size classes and with different eta and b values (based on the
# relationship between sinking speed and particle size)
fig = plt.figure(1, figsize=(3.5, 3.5))
ax = fig.add_axes([0.18, 0.18, width, height])
plt.boxplot([Flux_extended_eta_b_filtered_depthf,Flux_Clements_extended_filtered_depthf,Flux_extended_eta_b_filtered_depthf_231,POC_in_domain])
plt.ylim(0,70)
plt.ylabel('POC Flux (mgC/m$^2/d$)', fontsize=fs)
plt.title('Argo Flux with small size classes, eta=0.62,b=66,\n between %d-%02d-%02d and %d-%02d-%02d' % (day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day), fontsize=8)
plt.xticks([1,2,3,4],['BGC Argo\n6903095 ','BGC Argo\n6903095\n(Clements\net al.)','Transect\nMSM60','Sediment traps\n(Mouw et al.)'], fontsize=fs)
ax.text(-0.15, 1.125, 'b', transform=ax.transAxes, fontsize=12, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
plt.savefig('../Plots/an32/04extended_eta_b_POCFlux_Argo_vs_literature_%d%02d%02dto%d%02d%02d_an32.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day) ,dpi=200)
plt.close()


from scipy.stats import ttest_ind, mannwhitneyu,shapiro

mannwhitneyu(POC_in_domain, Flux_filtered_depthf)
mannwhitneyu(POC_in_domain, Flux_extended_filtered_depthf)
mannwhitneyu(POC_in_domain, Flux_extended_filtered_depthf_231)
mannwhitneyu(POC_in_domain, Flux_filtered_depthf_231)
mannwhitneyu(Flux_filtered_depthf, Flux_filtered_depthf_231)
mannwhitneyu(Flux_extended_filtered_depthf, Flux_extended_filtered_depthf_231)

# two-sample t-test
# null hypothesis: the two groups have the same mean
# this test assumes the two groups have the same variance and that they are gaussianly distributed
ttest_ind(POC_in_domain, Flux_filtered_depthf)
ttest_ind(POC_in_domain, Flux_filtered_depthf)
ttest_ind(POC_in_domain, Flux_extended_filtered_depthf)
ttest_ind(POC_in_domain, Flux_extended_filtered_depthf_231)
ttest_ind(POC_in_domain, Flux_filtered_depthf_231)
ttest_ind(Flux_filtered_depthf, Flux_filtered_depthf_231)
ttest_ind(Flux_extended_filtered_depthf, Flux_extended_filtered_depthf_231)
# To test if the distribution is guassian
shapiro(POC_in_domain) #not gaussian
shapiro(Flux_filtered_depthf) #not gaussian
shapiro(Flux_filtered_depthf_231) #not gaussian
shapiro(Flux_extended_filtered_depthf) #gaussian
shapiro(Flux_extended_filtered_depthf_231) #not gaussian
