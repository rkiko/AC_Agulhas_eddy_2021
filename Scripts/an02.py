import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename='../Data/Ecopart_mip_map_flux_data.tsv'
data=pd.read_csv(filename, sep='\t', header=0)
data.columns = data.columns.str.replace(' ','_') # I remove spaces and [] symbols
data.columns = data.columns.str.replace('[','')
data.columns = data.columns.str.replace(']','')
RAWfilename=data.RAWfilename

#I select only the prophiles data, which contain 'ASC' in the filename, and I exclude the parkings
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a.split('-')[-1].split('_')[0] == 'ASC':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1

# I extract the data
lon=np.array(data.Longitude[sel_filename]);lat=np.array(data.Latitude[sel_filename])
Flux=np.array(data.Flux_mgC_m2[sel_filename]);Date_Time=np.array(data.Date_Time[sel_filename])
pressure=-np.array(data.Pressure_dbar[sel_filename])
MiP=np.array(data.MiP_abun[sel_filename])
MaP=np.array(data.MaP_abun[sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I filter the flux prophiles
Flux_filtered=np.array([]);pressure_Flux=np.array([]);Date_Num_Flux=np.array([])
MiP_filtered=np.array([]);pressure_MiP=np.array([]);Date_Num_MiP=np.array([])
MaP_filtered=np.array([]);pressure_MaP=np.array([]);Date_Num_MaP=np.array([])
list_dates=np.unique(Date_Num)
for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i];x=Date_Num[sel];y=pressure[sel]
    z=Flux[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        Flux_filtered = np.concatenate((Flux_filtered, z));Date_Num_Flux = np.concatenate((Date_Num_Flux, x2));pressure_Flux = np.concatenate((pressure_Flux, y2))
    z=MiP[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        MiP_filtered = np.concatenate((MiP_filtered, z));Date_Num_MiP = np.concatenate((Date_Num_MiP, x2));pressure_MiP = np.concatenate((pressure_MiP, y2))
    z=MaP[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        MaP_filtered = np.concatenate((MaP_filtered, z));Date_Num_MaP = np.concatenate((Date_Num_MaP, x2));pressure_MiP = np.concatenate((pressure_MaP, y2))

# I define the x and y arrays for the contourf plot
x_Flux = np.linspace(Date_Num_Flux.min(),Date_Num_Flux.max(),100)
y_Flux = np.linspace(pressure_Flux.min(),pressure_Flux.max(),50)
x_MiP = np.linspace(Date_Num_MiP.min(),Date_Num_MiP.max(),100)
y_MiP = np.linspace(pressure_MiP.min(),pressure_MiP.max(),50)
x_Flux_g,y_Flux_g=np.meshgrid(x_Flux,y_Flux)
x_MiP_g,y_MiP_g=np.meshgrid(x_MiP,y_MiP)
# I interpolate
Flux_interp = griddata((Date_Num_Flux,pressure_Flux), Flux_filtered, (x_Flux_g, y_Flux_g), method="nearest")
MiP_interp = griddata((Date_Num_MiP,pressure_MiP), MiP_filtered, (x_MiP_g, y_MiP_g), method="nearest")

################################################################################################################
####### I plot
################################################################################################################

########################################################
####### Flux
########################################################
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -600, pressure_Flux.max()
#Scatter plot without filter
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
maxFlux=30;Flux_plot=Flux
Flux_plot[Flux_plot>maxFlux]=maxFlux
plot1=plt.scatter(Date_Num,pressure,c=np.log(Flux_plot))#, vmin=1, vmax=np.log(maxFlux))
cbar = fig.colorbar(plot1)
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
#plt.clim(np.log(1),np.log(maxFlux))
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxFlux),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test02_Flux_noFilter_an02.pdf',dpi=200)
#plt.show()
plt.close()

#Scatter plot of the filtered data
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
Flux_plot=Flux_filtered
Flux_plot[Flux_plot>maxFlux]=maxFlux
plot1=plt.scatter(Date_Num_Flux,pressure_Flux,c=np.log(Flux_plot))
cbar = fig.colorbar(plot1)
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxFlux),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test02_Flux_yesFilter_an02.pdf',dpi=200)
#plt.show()
plt.close()

####### I do the interpolate contourf
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -600, pressure_Flux.max()
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
Flux_plot=Flux_interp
Flux_plot[Flux_plot>maxFlux]=maxFlux
ax_1 = plot2 = plt.contourf(x_Flux, y_Flux, np.log(Flux_plot))#, vmax=np.log(maxFlux))#, cmap=cmhot)
# draw colorbar
cbar = plt.colorbar(plot2)
#plt.clim(np.log(1),np.log(maxFlux))
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxFlux),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test03_Flux_an02.pdf',dpi=200)
plt.close()



####### I do the interpolate contourf up to the maximum depth
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -600, pressure_Flux.max()
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(pressure_Flux.min(), set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
Flux_plot=Flux_interp
Flux_plot[Flux_plot>maxFlux]=maxFlux
ax_1 = plot2 = plt.contourf(x_Flux, y_Flux, np.log(Flux_plot))#, vmax=np.log(maxFlux))#, cmap=cmhot)
# draw colorbar
cbar = plt.colorbar(plot2)
#plt.clim(np.log(1),np.log(maxFlux))
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxFlux),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test03B_Flux_an02.pdf',dpi=400)
#plt.show()
plt.close()

########################################################
####### MiP
########################################################

maxMiP=30
####### I do the interpolate contourf
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -600, pressure_MiP.max()
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
MiP_plot=MiP_interp
MiP_plot[MiP_plot>maxMiP]=maxMiP
ax_1 = plot2 = plt.contourf(x_MiP, y_MiP, np.log(MiP_plot))#, vmax=np.log(maxMiP))#, cmap=cmhot)
# draw colorbar
cbar = plt.colorbar(plot2)
#plt.clim(np.log(1),np.log(maxFlux))
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('MiP (mg $m^{-3}$)', fontsize=18)
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxMiP),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test04_MiP_an02.pdf',dpi=200)
plt.close()

####### I do the interpolate contourf down to the maximum depth
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = pressure_MiP.min(), pressure_MiP.max()
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
MiP_plot=MiP_interp
MiP_plot[MiP_plot>maxMiP]=maxMiP
ax_1 = plot2 = plt.contourf(x_MiP, y_MiP, np.log(MiP_plot))#, vmax=np.log(maxMiP))#, cmap=cmhot)
# draw colorbar
cbar = plt.colorbar(plot2)
#plt.clim(np.log(1),np.log(maxFlux))
cbar.ax.get_yticklabels()
cbar.ax.set_ylabel('MiP (mg $m^{-3}$)', fontsize=18)
#I set ticks of the colorbar
ncbarticks=4
cbarticks=np.linspace(1,np.log(maxMiP),ncbarticks)
cbarticklabels=[]
for i in cbarticks:
    cbarticklabels.append('%d' % np.round(np.e**i))
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarticklabels)
plt.ylabel('Pressure (dbar)', fontsize=18)
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
    xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=12)
# I add the grid
plt.grid(color='k', linestyle='dashed', linewidth=0.5)
plt.savefig('../Plots/test04B_MiP_an02.pdf',dpi=200)
plt.close()

