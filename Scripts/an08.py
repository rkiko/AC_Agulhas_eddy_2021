import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import seawater as sw
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename_ecopart='%s/GIT/AC_Agulhas_eddy_2021/Data/Ecopart_mip_map_flux_data.tsv' % home
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
lon=np.array(data['Longitude'][sel_filename])
lat=np.array(data['Latitude'][sel_filename])
Date_Time=np.array(data['Date_Time'][sel_filename])
pressure=np.array(data['Pressure [dbar]'][sel_filename])
Flux=np.array(data['Flux_mgC_m2'][sel_filename])
MiP_abund=np.array(data['MiP_abun'][sel_filename])
MaP_abund=np.array(data['MaP_abun'][sel_filename])
MiP_POC=np.array(data['Mip_POC_cont_mgC_m3'][sel_filename])
MaP_POC=np.array(data['Map_POC_cont_mgC_m3'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

#I tranform the pressure to depth
mask_depth=pressure!=99999 #I select only valid values
lat_tmp=lat[mask_depth]
pressure_tmp=pressure[mask_depth]
depth_tmp=sw.eos80.dpth(pressure_tmp, lat_tmp)
depth=np.ones(pressure.shape)*99999
depth[mask_depth]=depth_tmp

list_dates=np.unique(Date_Num)
#I start the loop on the different parameters I plot
parameter=Flux
ipar=1
parameter_shortname_list=['Flux','MiP_abund','MaP_abund','MiP_POC','MaP_POC']
parameter_ylabel_list=['Flux (mgC $m^{-2}$ $d^{-1}$)','MiP abundance (# L$^{-1}$)','MaP abundance (# L$^{-1}$)'
    ,'MiP (mgC $m^{-3}$)','MaP (mgC $m^{-3}$)']
max_parameter_list=np.array([6,65,0.6,1.15,0.012])
for ipar in range(0,parameter_ylabel_list.__len__()):
    if ipar==0: parameter=Flux.copy()
    elif ipar==1:   parameter=MiP_abund.copy()
    elif ipar == 2: parameter=MaP_abund.copy()
    elif ipar == 3: parameter=MiP_POC.copy()
    elif ipar == 4: parameter=MaP_POC.copy()

    # I filter the flux prophiles
    parameter_filtered=np.array([]);depth_filtered=np.array([]);Date_Num_filtered=np.array([])
    i=0
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
        z=parameter[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2)>0:
            z=savgol_filter(z,5,1)
            parameter_filtered = np.concatenate((parameter_filtered, z))
            Date_Num_filtered = np.concatenate((Date_Num_filtered, x2))
            depth_filtered = np.concatenate((depth_filtered, y2))

    # I define the x and y arrays for the contourf plot
    x_filtered = np.linspace(Date_Num_filtered.min(),Date_Num_filtered.max(),100)
    y_filtered = np.linspace(depth_filtered.min(),depth_filtered.max(),100)
    x_filtered_g,y_filtered_g=np.meshgrid(x_filtered,y_filtered)
    # I interpolate
    parameter_interp = griddata((Date_Num_filtered,depth_filtered), parameter_filtered, (x_filtered_g, y_filtered_g), method="nearest")

    ################################################################################################################
    ####### I plot
    ################################################################################################################
    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),600
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot<0]=0
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    # draw colorbar
    cbar = plt.colorbar(plot2)
    cbar.ax.get_yticklabels()
    cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
    #I set ticks of the colorbar
    # ncbarticks=6
    # cbarticks=np.linspace(parameter_plot.min(),np.log(max_parameter_list[ipar]),ncbarticks)
    # cbarticklabels=[]
    # for i in cbarticks:
    #     cbarticklabels.append('%d' % np.round(np.e**i))
    # cbar.set_ticks(cbarticks)
    # cbar.set_ticklabels(cbarticklabels)
    plt.ylabel('Depth (m)', fontsize=18)
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
    plt.savefig('../Plots/an08/TimeSeries%02dA%s_an08.pdf' % (ipar,parameter_shortname_list[ipar]),dpi=200)
    plt.close()



    ####### I do the interpolate contourf up to the maximum depth
    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = depth_filtered.min(),1000
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    parameter_plot=parameter_interp
    parameter_plot[parameter_plot>max_parameter_list[ipar]]=max_parameter_list[ipar]
    ax_1 = plot2 = plt.contourf(x_filtered, y_filtered, parameter_plot)
    plt.gca().invert_yaxis()
    # draw colorbar
    cbar = plt.colorbar(plot2)
    cbar.ax.get_yticklabels()
    cbar.ax.set_ylabel(parameter_ylabel_list[ipar], fontsize=18)
    #I set ticks of the colorbar
    # ncbarticks=6
    # cbarticks=np.linspace(parameter_plot.min(),np.log(max_parameter_list[ipar]),ncbarticks)
    # cbarticklabels=[]
    # for i in cbarticks:
    #     cbarticklabels.append('%d' % np.round(np.e**i))
    # cbar.set_ticks(cbarticks)
    # cbar.set_ticklabels(cbarticklabels)
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
    plt.savefig('../Plots/an08/TimeSeries%02dB%s_an08.pdf' % (ipar,parameter_shortname_list[ipar]),dpi=200)
    plt.close()


