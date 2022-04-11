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


filename='../Data/Ecopart_diagnostics_data_356.tsv'
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
Date_Time=np.array(data.Date_Time[sel_filename]);
depth=np.array(data['Depth_m'][sel_filename])

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:depth.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

threshold_min_list=np.array([3])
threshold_max_list=np.array([8])
Date_Num_limit_min=1626742080+86400*np.array([0])
Date_Num_limit_max=1630358605+86400*np.array([-10])
nthresholds=5 #how many thresholds I consider for each size class
settling_vel_mean=np.array([]);settling_vel_std=np.array([]);sizeclass_list=np.array([])
# I define a list containing the index of the columns corresponding with the size classes for which I calculate the sinking
# velocity
list_indexes=[7]
idx=0
for idx in range(0, len(list_indexes)):
    iNP = list_indexes[idx]
    #########Begin text to indent
    threshold_value_list=np.linspace(threshold_min_list[idx],threshold_max_list[idx],nthresholds)
    Date_Num_limit=np.array([Date_Num_limit_min[idx],Date_Num_limit_max[idx]])
    NP_abun=data.values[sel_filename,iNP]
    NP_abun=NP_abun.astype('float')
    # NP_sizeclass=data.axes[1][iNP]
    # NP_sizeclass=NP_sizeclass.split('_')[-2]
    # sizeclass_list=np.append(sizeclass_list,(float(NP_sizeclass.split('-')[0])+float(NP_sizeclass.split('-')[1]))*0.5)
    settling_vel=np.array([])
    threshold_value=threshold_value_list[0]
    for threshold_value in threshold_value_list:
        # I filter the prophiles
        NP_filtered=np.array([]);depth_NP=np.array([]);Date_Num_NP=np.array([])
        list_dates=np.unique(Date_Num)
        depth_threshold=np.array([]);Date_Num_threshold=np.array([])
        for i in range(0,list_dates.size):
            sel=Date_Num==list_dates[i];x=Date_Num[sel];y=depth[sel]
            z=NP_abun[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
            if sum(sel2)>0:
                z=savgol_filter(z,5,1)
                NP_filtered = np.concatenate((NP_filtered, z));Date_Num_NP = np.concatenate((Date_Num_NP, x2));depth_NP = np.concatenate((depth_NP, y2))
                ct = 0
                for iz in range(0,z.size-1):
                    if (Date_Num_limit[0]<=list_dates[i]<=Date_Num_limit[1])&(z[iz+1]<threshold_value<z[iz])&(list_dates[i]!=1629078553): #I exclude an outlier
                        ct=ct+1
                        dist=(threshold_value-z[iz])/(z[iz+1]-z[iz])
                        depth_tmp=y[iz]-abs((y[iz+1]-y[iz])*dist)
                        #depth_threshold = np.append(depth_threshold, depth_tmp)
                        #Date_Num_threshold = np.append(Date_Num_threshold, list_dates[i])
                if ct==1:
                    depth_threshold = np.append(depth_threshold, depth_tmp)
                    Date_Num_threshold = np.append(Date_Num_threshold, list_dates[i])

        (interpol,_,_,signif,signif_label)=lin_fit(Date_Num_threshold/86400,depth_threshold)
        settling_vel=np.append(settling_vel,interpol.slope)
        print('Flux, Threshold %0.1f #/L, settling velocity: %0.2f m/d' % (threshold_value,interpol.slope))
        # I define the x and y arrays for the contourf plot
        x_NP = np.linspace(Date_Num_NP.min(), Date_Num_NP.max(), 100)
        y_NP = np.linspace(depth_NP.min(), depth_NP.max(), 50)
        x_NP_g, y_NP_g = np.meshgrid(x_NP, y_NP)
        # I interpolate
        NP_interp = griddata((Date_Num_NP, depth_NP), NP_filtered, (x_NP_g, y_NP_g), method="nearest")

        width, height = 0.8, 0.7
        set_ylim_lower, set_ylim_upper = 0, 600
        fig = plt.figure(1, figsize=(12, 8))
        ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper),
                          xlim=(Date_Num.min(), Date_Num.max()))
        NP_plot = NP_interp.copy()
        NP_plot[NP_plot>16]=16
        ax_1 = plot2 = plt.contourf(x_NP, y_NP, NP_plot)  # , vmax=np.log(maxNP))#, cmap=cmhot)
        # draw colorbar
        cbar = plt.colorbar(plot2)
        plot3 = plt.scatter(Date_Num_threshold, depth_threshold, c='black')
        plot4= plt.plot(np.linspace(Date_Num_limit[0],Date_Num_limit[1],20),np.linspace(Date_Num_limit[0]*interpol.slope/86400+interpol.intercept,Date_Num_limit[1]*interpol.slope/86400+interpol.intercept,20),c='black')
        cbar.ax.set_ylabel('Particle concentration (#/L)', fontsize=18)
        plt.gca().invert_yaxis()
        plt.ylabel('Depth (m)', fontsize=18)
        plt.title('Threshold: %0.1f mgC/m2/d, Settling vel: %0.2f m/d, fit signif: %s' % (threshold_value,interpol.slope,signif_label), fontsize=18)
        nxticks = 10
        xticks = np.linspace(Date_Num.min(), Date_Num.max(), nxticks)
        xticklabels = []
        for i in xticks:
            xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        plt.xticks(rotation=90, fontsize=12)
        ax.text(-0.07, 1.02, 'a', transform=ax.transAxes, fontsize=30, fontweight='bold',va='top', ha='right')  # ,fontfamily='helvetica'
        # I add the grid
        plt.grid(color='k', linestyle='dashed', linewidth=0.5)
        plt.savefig('../Plots/an51/Flux_slope_thr%0.1f_an51.pdf' % (threshold_value), dpi=200)
        #plt.show()
        #input("Press Enter to continue...")
        plt.close()
        #exit()
        #########End text to indent

    settling_vel_mean=np.append(settling_vel_mean,np.mean(settling_vel))
    settling_vel_std=np.append(settling_vel_std,np.std(settling_vel))

    print('Flux, mean settling velocity: %0.2f±%0.2f m/d' % (settling_vel_mean, settling_vel_std))


