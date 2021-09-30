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

filename='../Data/Ecopart_processed_data.tsv'
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
Date_Time=np.array(data.Date_Time[sel_filename]);pressure=-np.array(data.Pressure_dbar[sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:pressure.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

iNP=25
for iNP in range(25,39):
    NP_abun=data.values[sel_filename,iNP]
    NP_abun=NP_abun.astype('float')
    NP_sizeclass=data.axes[1][iNP]
    NP_sizeclass=NP_sizeclass.split('_')[-2]

    # I filter the prophiles
    NP_filtered=np.array([]);pressure_NP=np.array([]);Date_Num_NP=np.array([])
    list_dates=np.unique(Date_Num)
    for i in range(0,list_dates.size):
        sel=Date_Num==list_dates[i];x=Date_Num[sel];y=pressure[sel]
        z=NP_abun[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
        if sum(sel2)>0:
            z=savgol_filter(z,5,1)
            NP_filtered = np.concatenate((NP_filtered, z));Date_Num_NP = np.concatenate((Date_Num_NP, x2));pressure_NP = np.concatenate((pressure_NP, y2))

    # I define the x and y arrays for the contourf plot
    x_NP = np.linspace(Date_Num_NP.min(),Date_Num_NP.max(),100)
    y_NP = np.linspace(pressure_NP.min(),pressure_NP.max(),50)
    x_NP_g,y_NP_g=np.meshgrid(x_NP,y_NP)
    # I interpolate
    NP_interp = griddata((Date_Num_NP,pressure_NP), NP_filtered, (x_NP_g, y_NP_g), method="nearest")

    ################################################################################################################
    ####### I plot
    ################################################################################################################

    ########################################################
    ####### NP
    ########################################################

    #maxNP=np.round(NP_filtered.max()/4)
    ####### I do the interpolate contourf
    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = -600, pressure_NP.max()
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    NP_plot=NP_interp
    #NP_plot[NP_plot>maxNP]=maxNP
    ax_1 = plot2 = plt.contourf(x_NP, y_NP, np.log(NP_plot))#, vmax=np.log(maxNP))#, cmap=cmhot)
    # draw colorbar
    cbar = plt.colorbar(plot2)
    #plt.clim(np.log(1),np.log(maxNP))
    #a=cbar.ax.get_yticklabels()[-1];a=float(cbar.ax.get_yticklabels()[-1].get_text())
    cbar.ax.set_ylabel('NP abundance (#/L)', fontsize=18)
    #I set ticks of the colorbar
    ncbarticks=6
    cbarticks=np.linspace(cbar.ax.get_yticklabels()[0].get_position()[1],cbar.ax.get_yticklabels()[-1].get_position()[1],ncbarticks) # cbar.ax.get_yticklabels()[0].get_position()[1] is the maximum tick value in the colorbar
    cbarticklabels=[]
    for i in cbarticks:
        if iNP <= 27:
            cbarticklabels.append('%d' % np.round(np.e ** i))
        elif 27 < iNP <= 29 :
            cbarticklabels.append('%0.1f' % np.e ** i)
        elif 29 < iNP <= 31:
            cbarticklabels.append('%0.2f' % np.e ** i)
        else:
            cbarticklabels.append('%0.3f' % np.e ** i)
    cbar.set_ticks(cbarticks)
    cbar.set_ticklabels(cbarticklabels)
    plt.ylabel('Pressure (dbar)', fontsize=18)
    plt.title('%smm' % NP_sizeclass, fontsize=18)
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
    NP_sizeclass_save=NP_sizeclass
    if NP_sizeclass=='0.0806-0.102':    NP_sizeclass_save='0.081-0.102'
    plt.savefig('../Plots/an03/test05_NP%smm_an03.pdf' % NP_sizeclass_save,dpi=200)
    plt.close()



    ####### I do the interpolate contourf up to the maximum depth
    width, height = 0.8, 0.7
    set_ylim_lower, set_ylim_upper = pressure_NP.min(), pressure_NP.max()
    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
    NP_plot=NP_interp
    #NP_plot[NP_plot>maxNP]=maxNP
    ax_1 = plot2 = plt.contourf(x_NP, y_NP, np.log(NP_plot))#, vmax=np.log(maxNP))#, cmap=cmhot)
    # draw colorbar
    cbar = plt.colorbar(plot2)
    #plt.clim(np.log(1),np.log(maxNP))
    #a=cbar.ax.get_yticklabels()[-1];a=float(cbar.ax.get_yticklabels()[-1].get_text())
    cbar.ax.set_ylabel('NP abundance (#/L)', fontsize=18)
    #I set ticks of the colorbar
    ncbarticks=6
    cbarticks=np.linspace(cbar.ax.get_yticklabels()[0].get_position()[1],cbar.ax.get_yticklabels()[-1].get_position()[1],ncbarticks) # cbar.ax.get_yticklabels()[0].get_position()[1] is the maximum tick value in the colorbar
    cbarticklabels=[]
    for i in cbarticks:
        if iNP <= 27:
            cbarticklabels.append('%d' % np.round(np.e ** i))
        elif 27 < iNP <= 29 :
            cbarticklabels.append('%0.1f' % np.e ** i)
        elif 29 < iNP <= 31:
            cbarticklabels.append('%0.2f' % np.e ** i)
        else:
            cbarticklabels.append('%0.3f' % np.e ** i)
    cbar.set_ticks(cbarticks)
    cbar.set_ticklabels(cbarticklabels)
    plt.ylabel('Pressure (dbar)', fontsize=18)
    plt.title('%smm' % NP_sizeclass, fontsize=18)
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
    NP_sizeclass_save=NP_sizeclass
    if NP_sizeclass=='0.0806-0.102':    NP_sizeclass_save='0.081-0.102'
    plt.savefig('../Plots/an03/test05B_NP%smm_an03.pdf' % NP_sizeclass_save,dpi=200)
    plt.close()

