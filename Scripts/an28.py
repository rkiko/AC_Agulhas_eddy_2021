import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from lin_fit import lin_fit

########################################################################################################################
#Time period that I plot
########################################################################################################################
day0=datetime(2021,4,13)        # starting date
dayf=datetime(2021,8,19)        # final date
day0_float = calendar.timegm(day0.timetuple())
dayf_float = calendar.timegm(dayf.timetuple())
ndays = (dayf - day0).days  # number of days
delta_bin=20 #thickness (in meters) of the bin I used to calculate the psd distribution
depth_f=1000 #final depth (in meters) at which I calculate the psd distribution

list_depths=np.r_[0:depth_f+delta_bin*0.5:delta_bin]
########################################################################################################################
#I process the data
########################################################################################################################
data = pd.read_csv("../Data/Ecopart_processed_data.tsv", sep='\t')

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
PSD_columns=data.columns[28:42]

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:depth.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# list_dates=np.sort(np.unique(Date_Num))
# list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
sel_dates=(Date_Num>=day0_float)&(Date_Num<=dayf_float)
depth=depth[sel_dates]

########################################################################################################################
#I initialise the PSD vs depth as a 2D matrix with a number of rows equal to the number of depth bins, and a number of
#columns equal to the number of size classes selected (between 0.1 and 2.6 mm, see above). I also initialise the size of
#the bins which I need to normalise the PSD distribution
########################################################################################################################
PSD_depth=np.zeros((list_depths.size-1,PSD_columns.size))
PSD_bin_start=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_end = np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_width=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_mean=np.squeeze(np.zeros((PSD_columns.size,1)))

########################################################################################################################
#I do the loop for each size class and each depth bin
########################################################################################################################

i_PSD = 0
for i_PSD in range(0,PSD_columns.size):
    PSD_columns_name = PSD_columns[i_PSD]
    PSD_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_bin_width[i_PSD] = PSD_bin_end[i_PSD] - PSD_bin_start[i_PSD]
    PSD_bin_mean[i_PSD] = (PSD_bin_end[i_PSD] + PSD_bin_start[i_PSD])*0.5
    PartConc = np.array(data[PSD_columns_name][sel_filename])
    PartConc = PartConc[sel_dates]

    i_depth=0
    for i_depth in range(0,list_depths.size-1):
        depth0_tmp = list_depths[i_depth]
        depthf_tmp = list_depths[i_depth + 1]
        sel_depth=(depth>=depth0_tmp)&(depth<depthf_tmp)
        PartConc_depth=PartConc[sel_depth]
        PSD_depth[i_depth,i_PSD]=np.nansum(PartConc_depth)/PSD_bin_width[i_PSD] #I divide by PSD_bin_width[i_PSD] in order to normalise

########################################################################################################################
#I plot, once for every depth bin
########################################################################################################################

i_depth=0
for i_depth in range(0,list_depths.size-1):
    depth0_tmp = list_depths[i_depth]
    depthf_tmp = list_depths[i_depth + 1]
    depth_tmp = (depthf_tmp + depth0_tmp) * 0.5
    PSD_depth_tmp=PSD_depth[i_depth,:]
    PSD_depth_tmp=PSD_depth_tmp
    x = np.log(PSD_bin_mean)
    y = np.log(PSD_depth_tmp)
    sel=~np.isinf(y);x=x[sel];y=y[sel]
    (interpol, slpe_ci, _, signif, signif_label) = lin_fit(x,y)

    width, height = 0.7, 0.68
    set_ylim_lower, set_ylim_upper = 0.1, 10**12
    set_xlim_lower, set_xlim_upper = 0.001,2.7
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.2, 0.25, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(set_xlim_lower, set_xlim_upper))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(np.exp(np.linspace(np.log(set_xlim_lower), set_xlim_upper, 20)),
             np.exp(np.linspace(np.log(set_xlim_lower) * interpol.slope + interpol.intercept,
                                set_xlim_upper * interpol.slope + interpol.intercept, 20)), c='red',label='Fit')
    plt.plot(PSD_bin_mean, PSD_depth_tmp,c='blue',label='data')
    plt.xlabel('Size (mm)',fontsize=10)
    plt.ylabel('Normalised Abundance (#/L)',fontsize=10)
    plt.legend(fontsize=10)
    plt.title('Depth: %d m; R$^2$=%0.2f, Signif: %s' % (depth_tmp,interpol.rvalue**2,signif_label),fontsize=10)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an28/PSD_%d%02d%02dto%d%02d%02d_depth%dm_an28.pdf' % (day0.year,day0.month,day0.day,dayf.year,dayf.month,dayf.day,depth_tmp),dpi=200)
    plt.close()

