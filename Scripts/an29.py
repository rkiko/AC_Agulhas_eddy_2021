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
# day0=datetime(2021, 4,13)        # starting date
# dayf=datetime(2021,10,21)        # final date
# day0_float = calendar.timegm(day0.timetuple())
# dayf_float = calendar.timegm(dayf.timetuple())
# ndays = (dayf - day0).days  # number of days

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
ic0 = 28                                          # Starting index of the columns containing the size classes for which UVP works, and that I use for the fit
icf = 42                                          # End index of the columns containing the size classes for which UVP works, and that I use for the fit
ice0 = 12                                         # Starting index of the columns containing the size classes for which UVP does not works, and for which I extrapolate the abundance
icef = ic0                                        # End index of the columns containing the size classes for which UVP does not works, and for which I extrapolate the abundance
PSD_columns=data.columns[ic0:icf]                 # Column labels containing the size classes for which UVP works, and that I use for the fit
PSD_columns_extrapolation=data.columns[ice0:icef] # Column labels containing the size classes for which UVP does not works, and for which I extrapolate the abundance

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:depth.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

list_dates=np.sort(np.unique(Date_Num))
# list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
# sel_dates=(Date_Num>=day0_float)&(Date_Num<=dayf_float)
# depth=depth[sel_dates]

########################################################################################################################
#I copy the concentration for the size class 0.0806-0.102 mm in the last column because I want to keep this values for
#the extraction of the sinking speed
########################################################################################################################

data['Non-extrapol. Part conc frac [#/l] (ESD: 0.0806-0.102 mm)'] = data.values[:,27].copy()

########################################################################################################################
#I calculate the size of the bins which I need to normalise the PSD distribution
########################################################################################################################

PSD_bin_start=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_end = np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_width=np.squeeze(np.zeros((PSD_columns.size,1)))
PSD_bin_mean=np.squeeze(np.zeros((PSD_columns.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns.size):
    PSD_columns_name = PSD_columns[i_PSD]
    PSD_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_bin_width[i_PSD] = PSD_bin_end[i_PSD] - PSD_bin_start[i_PSD]
    PSD_bin_mean[i_PSD] = (PSD_bin_end[i_PSD] + PSD_bin_start[i_PSD])*0.5

########################################################################################################################
#I extract the width and mean values for the smaller size classes
########################################################################################################################

PSD_extrapolated_bin_start=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_end = np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_width=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))
PSD_extrapolated_bin_mean=np.squeeze(np.zeros((PSD_columns_extrapolation.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns_extrapolation.size):
    PSD_columns_name = PSD_columns_extrapolation[i_PSD]
    PSD_extrapolated_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_extrapolated_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_extrapolated_bin_width[i_PSD] = PSD_extrapolated_bin_end[i_PSD] - PSD_extrapolated_bin_start[i_PSD]
    PSD_extrapolated_bin_mean[i_PSD] = (PSD_extrapolated_bin_end[i_PSD] + PSD_extrapolated_bin_start[i_PSD])*0.5

PSD_extrapolated_bin_mean_log=np.log(PSD_extrapolated_bin_mean)

########################################################################################################################
#I do the loop for each profile, each depth, and each size class
########################################################################################################################
PSD_extrapolated=data.values[:,ice0:icef].copy()

ct_ok=0;ct_notOk=0
i_profile = 0
for i_profile in range (0,list_dates.size):
    sel_dates = Date_Num == list_dates[i_profile]
    list_depths = depth[sel_dates]
    i_depth = 1
    for i_depth in range(0, list_depths.size):
        PSD_depth = np.squeeze(np.zeros((PSD_columns.size,1)))
        i_PSD = 0
        for i_PSD in range(0, PSD_columns.size):
            PSD_columns_name = PSD_columns[i_PSD]
            PartConc = np.array(data[PSD_columns_name][sel_filename])  # I select a given size class data and I exclude the parkings
            PartConc = PartConc[sel_dates]                             # I select the data for the specific profile (i.e., date)
            PSD_depth[i_PSD] = PartConc[i_depth]                         # I select the data for the specific depth

        PSD_depth = PSD_depth/PSD_bin_width             # I divide by PSD_bin_width in order to normalise
        sel = (~np.isnan(PSD_depth)) & (PSD_depth != 0) # I exlcude nan and zero values
        PSD_depth=PSD_depth[sel]
        #I extrapolate the concentration only if the size distribution has at least 3 values
        if PSD_depth.size>2:
            PSD_bin_mean_tmp = PSD_bin_mean[sel]
            x = np.log(PSD_bin_mean_tmp)
            y = np.log(PSD_depth)
            sel=~np.isinf(y);x=x[sel];y=y[sel]              # I exclude inf values
            (interpol, slpe_ci, _, signif, signif_label) = lin_fit(x,y)
            #I extrapolate the concentration only if the fit is significant or the pvalue is elss than 0.05
            if (signif>0)|(interpol.pvalue<=0.05):
                ct_ok = ct_ok +1
                PSD_extrapolated_depth_tmp = np.exp(PSD_extrapolated_bin_mean_log * interpol.slope + interpol.intercept)
                PSD_extrapolated[ np.squeeze(np.where(sel_filename))[sel_dates][i_depth] ,: ] = PSD_extrapolated_depth_tmp*PSD_extrapolated_bin_width
            else:
                ct_notOk = ct_notOk +1
        else:
            ct_notOk = ct_notOk + 1

PSD_extrapolated = PSD_extrapolated.astype('float32')

for i in range(0, PSD_columns_extrapolation.__len__()):
    data[PSD_columns_extrapolation[i]] = PSD_extrapolated[:, i]




########################################################################################################################
########################################################################################################################
########################################################################################################################
#I convert each size class concentration in POC and flux
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
#I import the functions necessary
########################################################################################################################

from paruvpy import poc_flux_func_1_class
from paruvpy import poc_cont_func_1_class

########################################################################################################################
#Time period during which I want to calculate the POC and flux and depth bins to average it
########################################################################################################################
day0=datetime(2021, 4,13)        # starting date
dayf=datetime(2021,8,19)        # final date
day0_float = calendar.timegm(day0.timetuple())
dayf_float = calendar.timegm(dayf.timetuple())
ndays = (dayf - day0).days  # number of days

delta_bin=20 #thickness (in meters) of the bin I used to calculate the psd distribution
depth_f=1000 #final depth (in meters) at which I calculate the psd distribution

list_depths=np.r_[0:depth_f+delta_bin*0.5:delta_bin]

# list_dates=np.sort(np.unique(Date_Num))
list_dates=list_dates[(list_dates>=day0_float)&(list_dates<=dayf_float)]
sel_dates=(Date_Num>=day0_float)&(Date_Num<=dayf_float)
depth=depth[sel_dates]

########################################################################################################################
#I extract the width and mean values for the all the size classes
########################################################################################################################
PSD_columns_all=data.columns[ice0:icf]

PSD_all_bin_start=np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_end = np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_width=np.squeeze(np.zeros((PSD_columns_all.size,1)))
PSD_all_bin_mean=np.squeeze(np.zeros((PSD_columns_all.size,1)))

i_PSD = 0
for i_PSD in range(0,PSD_columns_all.size):
    PSD_columns_name = PSD_columns_all[i_PSD]
    PSD_all_bin_start[i_PSD]=float(PSD_columns_name.split()[-2].split('-')[0])
    PSD_all_bin_end[i_PSD] = float(PSD_columns_name.split()[-2].split('-')[1])
    PSD_all_bin_width[i_PSD] = PSD_all_bin_end[i_PSD] - PSD_all_bin_start[i_PSD]
    PSD_all_bin_mean[i_PSD] = (PSD_all_bin_end[i_PSD] + PSD_all_bin_start[i_PSD])*0.5

########################################################################################################################
# I calculate the POC and Flux and plot them. For the POC, I add the bbp POC for comparison
########################################################################################################################
bbp=np.array(data['bbp POC [mgC/m3]'][sel_filename][sel_dates])

i_depth=0
for i_depth in range(0,list_depths.size-1):
    depth0_tmp = list_depths[i_depth]
    depthf_tmp = list_depths[i_depth + 1]
    depth_tmp = (depthf_tmp + depth0_tmp) * 0.5
    PartConc = data[PSD_columns_all][sel_filename][sel_dates]
    sel_depth = (depth >= depth0_tmp) & (depth < depthf_tmp)
    PartConc_depth = PartConc[sel_depth]
    PartConc_depth = np.mean(PartConc_depth,0)
    bbp_depth = np.mean(bbp [sel_depth])
    bbp_depth_std = np.std(bbp [sel_depth])

    Flux=np.squeeze(np.zeros((1,PSD_columns_all.size)))
    POC=np.squeeze(np.zeros((1,PSD_columns_all.size)))
    i_PSD = 0
    for i_PSD in range(0,PSD_columns_all.size):
        PSD_columns_name = PSD_columns_all[i_PSD]
        Flux[i_PSD]= PartConc_depth[i_PSD] * poc_flux_func_1_class(PSD_columns_name)
        POC[i_PSD] = PartConc_depth[i_PSD] * poc_cont_func_1_class(PSD_columns_name)

    width, height = 0.75, 0.75
    set_ylim_lower, set_ylim_upper = 0.5*10**-3, 5*10 ** 3
    set_xlim_lower, set_xlim_upper = 0.001, 2.7
    ##############
    #POC plot
    ##############
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.2, 0.15, width, height], xlim=(set_xlim_lower, set_xlim_upper), ylim=(set_ylim_lower, set_ylim_upper))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(PSD_all_bin_mean, POC, c='blue', label='UVP')
    plt.scatter(0.016, bbp_depth, label='bbp',c='red')
    plt.errorbar(0.016, bbp_depth, xerr=0.014,yerr=bbp_depth_std, capsize=5,c='red')
    plt.xlabel('Size (mm)', fontsize=10)
    plt.ylabel('POC (mgC/m$^3$)', fontsize=8)
    plt.legend(fontsize=10)
    plt.title('Depth: %d m' % depth_tmp, fontsize=10)#; R$^2$=%0.2f, Signif: %s' % (depth_tmp, interpol.rvalue ** 2, signif_label)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an29/01_POC_vs_Size_%d%02d%02dto%d%02d%02d_depth%dm_an29.pdf' % (
    day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day, depth_tmp), dpi=200)
    plt.close()

    ##############
    #Flux plot
    ##############
    set_ylim_lower, set_ylim_upper = 0.5 * 10 ** -1, 5 * 10 ** 3
    fig = plt.figure(1, figsize=(3.5, 3.5))
    ax = fig.add_axes([0.2, 0.15, width, height], xlim=(set_xlim_lower, set_xlim_upper), ylim=(set_ylim_lower, set_ylim_upper))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(PSD_all_bin_mean, Flux, c='blue', label='data')
    plt.xlabel('Size (mm)', fontsize=10)
    plt.ylabel('Flux (mgC/m$^2$/d)', fontsize=8)
    # plt.legend(fontsize=10)
    plt.title('Depth: %d m' % depth_tmp, fontsize=10)#; R$^2$=%0.2f, Signif: %s' % (depth_tmp, interpol.rvalue ** 2, signif_label)
    plt.grid(color='k', linestyle='dashed', linewidth=0.5)
    plt.savefig('../Plots/an29/02_Flux_vs_Size_%d%02d%02dto%d%02d%02d_depth%dm_an29.pdf' % (
    day0.year, day0.month, day0.day, dayf.year, dayf.month, dayf.day, depth_tmp), dpi=200)
    plt.close()

