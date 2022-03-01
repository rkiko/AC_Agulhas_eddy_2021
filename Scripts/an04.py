import calendar
import numpy as np
import os
import pandas as pd
from datetime import date,datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import scipy.stats
import sys
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts/" % home)
from lin_fit import lin_fit


filename='../Data/Ecopart_processed_data_356.tsv'
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

threshold_min_list=np.array([340,12,7,3,1.,0.4])
threshold_max_list=np.array([400,21,13,6,2.5,.8])
Date_Num_limit_min=1626742080+86400*np.array([0,0,0,0,0,0])
Date_Num_limit_max=1630358605+86400*np.array([0,0,0,0,-14,-14])
nthresholds=5 #how many thresholds I consider for each size class
settling_vel_mean=np.array([]);settling_vel_std=np.array([]);sizeclass_list=np.array([])
iNP=30
for iNP in range(25,31):
#########Begin text to indent
    threshold_value_list=np.linspace(threshold_min_list[iNP-25],threshold_max_list[iNP-25],nthresholds)
    if iNP<=27: threshold_value_list=np.round(threshold_value_list)
    threshold_value=threshold_value_list[3]
    Date_Num_limit=np.array([Date_Num_limit_min[iNP-25],Date_Num_limit_max[iNP-25]])
    NP_abun=data.values[sel_filename,iNP]
    NP_abun=NP_abun.astype('float')
    NP_sizeclass=data.axes[1][iNP]
    NP_sizeclass=NP_sizeclass.split('_')[-2]
    sizeclass_list=np.append(sizeclass_list,(float(NP_sizeclass.split('-')[0])+float(NP_sizeclass.split('-')[1]))*0.5)
    settling_vel=np.array([])
    threshold_value=threshold_value_list[0]
    for threshold_value in threshold_value_list:
        # I filter the prophiles
        NP_filtered=np.array([]);pressure_NP=np.array([]);Date_Num_NP=np.array([])
        list_dates=np.unique(Date_Num)
        pressure_threshold=np.array([]);Date_Num_threshold=np.array([])
        for i in range(0,list_dates.size):
            sel=Date_Num==list_dates[i];x=Date_Num[sel];y=pressure[sel]
            z=NP_abun[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
            if sum(sel2)>0:
                z=savgol_filter(z,5,1)
                NP_filtered = np.concatenate((NP_filtered, z));Date_Num_NP = np.concatenate((Date_Num_NP, x2));pressure_NP = np.concatenate((pressure_NP, y2))
                ct = 0
                for iz in range(0,z.size-1):
                    if np.logical_and(Date_Num_limit[0]<=list_dates[i]<=Date_Num_limit[1],z[iz+1]<threshold_value<z[iz]):
                        ct=ct+1
                        dist=(threshold_value-z[iz])/(z[iz+1]-z[iz])
                        pressure_tmp=y[iz]-abs((y[iz+1]-y[iz])*dist)
                        #pressure_threshold = np.append(pressure_threshold, pressure_tmp)
                        #Date_Num_threshold = np.append(Date_Num_threshold, list_dates[i])
                if ct==1:
                    pressure_threshold = np.append(pressure_threshold, pressure_tmp)
                    Date_Num_threshold = np.append(Date_Num_threshold, list_dates[i])

        (interpol,_,_,signif,signif_label)=lin_fit(Date_Num_threshold/86400,pressure_threshold)
        settling_vel=np.append(settling_vel,interpol.slope)
        print('%s mm, Threshold %0.1f #/L, settling velocity: %0.2f m/d' % (NP_sizeclass,threshold_value,interpol.slope))
        # I define the x and y arrays for the contourf plot
        x_NP = np.linspace(Date_Num_NP.min(), Date_Num_NP.max(), 100)
        y_NP = np.linspace(pressure_NP.min(), pressure_NP.max(), 50)
        x_NP_g, y_NP_g = np.meshgrid(x_NP, y_NP)
        # I interpolate
        NP_interp = griddata((Date_Num_NP, pressure_NP), NP_filtered, (x_NP_g, y_NP_g), method="nearest")

        width, height = 0.8, 0.7
        set_ylim_lower, set_ylim_upper = -600, pressure_NP.max()
        fig = plt.figure(1, figsize=(12, 8))
        ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper),
                          xlim=(Date_Num.min(), Date_Num.max()))
        NP_plot = NP_interp
        # NP_plot[NP_plot>maxNP]=maxNP
        ax_1 = plot2 = plt.contourf(x_NP, y_NP, NP_plot)  # , vmax=np.log(maxNP))#, cmap=cmhot)
        # draw colorbar
        cbar = plt.colorbar(plot2)
        plot3 = plt.scatter(Date_Num_threshold, pressure_threshold, c='black')
        plot4= plt.plot(np.linspace(Date_Num_limit[0],Date_Num_limit[1],20),np.linspace(Date_Num_limit[0]*interpol.slope/86400+interpol.intercept,Date_Num_limit[1]*interpol.slope/86400+interpol.intercept,20),c='black')
        cbar.ax.set_ylabel('NP abundance (#/L)', fontsize=18)
        plt.ylabel('Pressure (dbar)', fontsize=18)
        plt.title('%smm, thr: %0.1f #/L, Settling vel: %0.2f m/d, fit signif: %s' % (NP_sizeclass,threshold_value,interpol.slope,signif_label), fontsize=18)
        nxticks = 10
        xticks = np.linspace(Date_Num.min(), Date_Num.max(), nxticks)
        xticklabels = []
        for i in xticks:
            xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        plt.xticks(rotation=90, fontsize=12)
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        NP_sizeclass_save = NP_sizeclass
        if NP_sizeclass == '0.0806-0.102':    NP_sizeclass_save = '0.081-0.102'
        plt.savefig('../Plots/an04/test_slope_NP%smm_thr%0.1f_an04.pdf' % (NP_sizeclass_save,threshold_value), dpi=200)
        #plt.show()
        #input("Press Enter to continue...")
        plt.close()
        #exit()
        #########End text to indent

    settling_vel_mean=np.append(settling_vel_mean,np.mean(settling_vel))
    settling_vel_std=np.append(settling_vel_std,np.std(settling_vel))

fig = plt.figure(1, figsize=(12, 8))
plt.plot(sizeclass_list,np.abs(settling_vel_mean),'o')
plt.plot(sizeclass_list,np.abs(settling_vel_mean))
plt.errorbar(sizeclass_list,np.abs(settling_vel_mean),yerr=settling_vel_std,capsize=5)
plt.grid(color='k', linestyle='dashed', linewidth=0.3)
plt.xlabel('Size class (mm)', fontsize=18)
plt.ylabel('Settling velocity (m/d)', fontsize=18)
plt.title('Relationship between size class and settling velocity', fontsize=18)
plt.xticks(fontsize=12),plt.yticks(fontsize=12)
plt.savefig('../Plots/an04/relationship_size_vs_settlingVelocity_an04.pdf', dpi=200)
plt.close()

