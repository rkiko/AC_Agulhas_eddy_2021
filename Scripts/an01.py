import calendar

import numpy as np
import os
import pandas as pd
from datetime import date,datetime
from calendar import timegm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from pathlib import Path
home = str(Path.home())
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

filename='../Data/Ecopart_mip_map_flux_data.tsv'
data=pd.read_csv(filename, sep='\t', header=0)#np.loadtxt(filename, delimiter=',', skiprows=1)
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

# I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:Flux.size]
for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    #datetime.utcfromtimestamp(Date_Num[i])

# I define the x and y arrays for the plot
x_date = np.linspace(Date_Num.min(),Date_Num.max(),100)
y_pressure = np.linspace(pressure.min(),pressure.max(),50)
x_date_g,y_pressure_g=np.meshgrid(x_date,y_pressure)

# I interpolate
Flux_interp = griddata((Date_Num,pressure), Flux, (x_date_g, y_pressure_g), method="linear")

####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = pressure.min(), pressure.max()
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_axes([0.12, 0.2, width, height], ylim=(set_ylim_lower, set_ylim_upper), xlim=(Date_Num.min(), Date_Num.max()))
ax_1 = plot2 = plt.contourf(x_date, y_pressure, Flux_interp)#, cmap=cmhot)
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Flux (mgC $m^{-2}$ $d^{-1}$)', fontsize=18)
plt.xlabel('Date', fontsize=18)
plt.ylabel('Pressure (dbar)', fontsize=18)
plt.title('Flux in AC Agulhas eddy', fontsize=18)
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
plt.savefig('../Plots/test01_Flux.png',dpi=200)

plt.show()
exit()



#! /usr/bin/env python

"""

in v3 new size limits used, leaving out the smallest

"""

import psycopg2 as pgsql
import pandas.io.sql as sql
import numpy as np
import pylab as pyl
import re
import os


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


from own_func import nlcmap #contains the nonlinear colormap definition
from own_func import sql_select_particle_abun_func, sql_select_poc_flux_func, sql_select_tricho_abun_func, sql_select_poc_cont_func
from module_uvp_plotting_v1 import *


from seawater import dpth

import argparse



################################################################################
######                           sql select                              #######
################################################################################


mydb = pgsql.connect(database="local_use", user="rkiko", password="rkiko")


################################################################################
######                      argument definition                          #######
################################################################################


cruise_choices = ["MSM022", "M106", "M119", "Cassiopee_1", 
                  "Cassiopee_2", "PS088b", "Tara_1", "Tara_2", "RB_p16n"]
                  
parameter_choices = ["ADCP", "Abun_a", "Abun_b", "Abun_c", "POC_hFlux", "POC_vFlux", "POC_cont", "Tricho", "Proto", "Zoo", "Cops"]

parser = argparse.ArgumentParser(description = 'Script to define equatorial snowfall figure')
parser.add_argument('-c1', '--cruise1', default = 'M106', type = str,\
 help = "List cruises to be plotted", choices = cruise_choices)
parser.add_argument('-c2', '--cruise2', default = 'M119', type = str,\
 help = "List cruises to be plotted", choices = cruise_choices)
parser.add_argument('-c3', '--cruise3', default = 'MSM022', type = str,\
 help = "List cruises to be plotted", choices = cruise_choices)
parser.add_argument('--measurement1', '-m1', default = 'ADCP', type = str,\
choices = parameter_choices, help = "define measurement to be plotted in leftmost plot")
parser.add_argument('--measurement2', '-m2', default = 'POC_cont', type = str,\
choices = parameter_choices, help = "define measurement to be plotted in leftmost plot")
parser.add_argument('--measurement3', '-m3', default = 'POC_hFlux', type = str,\
choices = parameter_choices, help = "define measurement to be plotted in leftmost plot")
parser.add_argument('--depth_select', '-ds', default = 5000, type = str,\
choices = ['5000', '1000'], help = "define scaling to use for 5000 m plots or 1000 m plots")


args = parser.parse_args()


depth_select = args.depth_select




################################################################################
######          define parameter, size resolution, depth tresholds,      #######
######           lon and lat, cruise, grid_method                        #######
################################################################################


#Define the instrument and Size resolution to be plotted:

instrument = "UVP" #can only be UVP

depth_tresh_lower = int(depth_select) #depth_threshold for selection
depth_tresh_upper = 0 #depth_threshold for selection


lat_n = 5.1
lat_s = -5.1

use_bottom_topo = 1 #choose if to use or not





################################################################################
######          define figure settings                                   #######
################################################################################




fig = plt.figure(1, figsize=(3.5,5.3))
#rect = l,b,w,h

#defining colormaps to be used:

levels1 = [0, 1, 2, 3, 6]

cmap_nonlin = nlcmap(pyl.cm.jet, levels1)
cmhot = plt.cm.get_cmap("seismic")

grid_method = "linear" #choose the gridding method here
#method to do the regridding, can be nearest (for nearest neighbour) or linear
#levels for non-linear plotting: [0, 1, 2, 3, 6, 10] are good for particle abundance along 23W

plt.rcParams.update({'font.size': 6})


set_ylim_lower = -int(depth_select)
set_ylim_upper = -10

xlimit_ = (-5,5)

width, height = 0.25, 0.4

plot_position_dict = {1:0.12,2:0.42,3: 0.72, 4:0.12, 5:0.42, 6:0.72}
plot_top_bottom_dict = {1:0.15, 2:0.15, 3:0.15, 4:0.58, 5:0.58, 6:0.58}

levels = 20.0 #levels for contour plotting


lat_grid = [-5,5]
depth_grid = [0,-5000]

################################################################################
#######     defining dictionaries for selection and plotting parameters  #######
################################################################################


measure_dict = {"POC_vFlux":[1,['a','e']]}

letter_dict = {6:'c', 
               5:'b', 
               4:'a', 
               3:'f', 
               2:'e', 
               1:'d', }


cruises_23W = ["'MSM022'", "'M106'", "'PS088b'", "'M119'", "'M130'"]


lon_dict = {  "Cassiopee_2":[164.5,165.5]
            , "Cassiopee_1":[157,158]
            , "Tara_2":[-160.,-150.]
            , "Tara_1":[-86., -83.]
            , "MSM022": [-23.5,-22.5]
            , "M106": [-23.5, -22.5]
            , "M119": [-23.5, -22.5]
            , "PS088b": [-23.5, -22.5]
            , "RB_p16n": [-160., -140.]
            , "M130": [-23.5, -22.5]}



                

if int(depth_select) == 5000:
    contour_def_dict = {  'ADCP':[-0.8,0.8]
    					, 'Abun_a':[0,100]
    					, 'Abun_b':[0,4.5]
    					, 'Abun_c':[0,20.]
    					, 'POC_hFlux':[-0.2,0.2]
    					, 'POC_vFlux':[0,8]
    					, 'POC_cont': [0,1]
    					, 'Tricho': [0,200]
    					, 'Zoo': [0,30]
    					, 'Cops': [0,30]
    					, 'Proto': [0,30]}
    					

elif int(depth_select) == 1000:					
    contour_def_dict = {  'ADCP':[-0.8,0.8]
    					, 'Abun_a':[0,100]
    					, 'Abun_b':[0,12.]
    					, 'Abun_c':[0,20.]
    					, 'POC_vFlux':[0,20.]
    					, 'POC_hFlux':[-1.,1.]
    					, 'POC_cont': [0,4]
    					, 'Tricho': [0,200]
    					, 'Zoo': [0,30]
    					, 'Cops': [0,30]
    					, 'Proto': [0,30]}
    					



################################################################################
####end of parameter definition, start of for loop for different measures#######
################################################################################


cruise_list = [args.cruise1, args.cruise2, args.cruise3]

cruise_list = ["M106", "M130", "M119", "PS088b", 
                  "MSM022", "RB_p16n"]


if cruise_list[0] == cruise_list[1]:
    cruise_list = [cruise_list[0]]


cruise_counter = 1 #set to 1 to select the bottom definition for plotting



for c in cruise_list:
    print(cruise_counter)
    
    ### general settings
    
    lon_w = lon_dict[c][0]
    lon_e = lon_dict[c][1]
    fig_level = plot_top_bottom_dict[cruise_counter]
    plot_position = plot_position_dict[cruise_counter]
    print(fig_level)
    
        
    if c.startswith("Cassiop"):
        cruise = "'Cassiopee'"
    elif  c.startswith("Tara"):
        cruise = "'Tara'"
    else:
        cruise = "'" + c + "'"
    
    ### topography
    
    
    if cruise in cruises_23W:
    
        df_topo = sql.read_sql(""" 
        select topo.PositionLat_dec, topo.Depth
        FROM topo_23w as topo
        where topo.PositionLat_dec BETWEEN -5 and 5
        ORDER by PositionLat_dec
        """, mydb)
        
        bottom_topo_lat =  df_topo['positionlat_dec'].tolist()
        bottom_topo_depth = (df_topo['depth']*-1).tolist()
    
    
    else:
        from own_func import transect_topo
        lon_topo = [lon_w, lon_e]
        topo = transect_topo(-5.,5., np.mean(lon_topo), np.mean(lon_topo))
        bottom_topo_depth, bottom_topo_lat = topo[0], topo[1]
    
    
    
    ### subplotting
    
    
    for item in measure_dict:
        measure = item
        print(measure)
        
        
        ### subplot defintions
        
        
        title = title_dict[item]
    
        
        if 'ADCP' in measure or 'hFlux' in measure:
            cmap_ = cmhot
        else:
            cmap_ = cmap_nonlin
        
        ### contour definitions
        
        min_contour_level = contour_def_dict[item][0]
        max_contour_level = contour_def_dict[item][1]
        distance_levels = max_contour_level/levels
        contour_levels = np.arange(min_contour_level,max_contour_level,distance_levels)
        
        
        
        if measure != "ADCP":
            measure_sql = sql_select_dict[measure][0]
            sql_select = sql_select_dict[measure][1]
        
        
        if measure in ['ADCP', 'POC_hFlux']:
            
            if c in ["M119", "RB_p16n", "PS088b", "Cassiopee_1", "Cassiopee_2"]:
                
                df_lADCP = sql.read_sql(""" 
            
                select lADCP.PositionLat_dec, lADCP.Depth, lADCP.u
                FROM lADCP_profile_data as lADCP
                WHERE lADCP.Cruise= %s AND lADCP.Depth BETWEEN 0 AND 6000 
                AND lADCP.PositionLat_dec BETWEEN %s AND %s
                AND lADCP.PositionLon_dec BETWEEN %s AND %s
                
                """%(cruise, lat_s, lat_n, lon_dict[c][0], lon_dict[c][1]), mydb)
                
                if cruise == "'RB_p16n'":
                    df_lADCP_min = sql.read_sql(""" 
                
                    select lADCP.PositionLat_dec as lat, -1*max(lADCP.Depth) as max_depth
                    FROM lADCP_profile_data as lADCP
                    WHERE lADCP.Cruise= %s AND lADCP.Depth BETWEEN 0 AND 6000 
                    AND lADCP.PositionLat_dec BETWEEN %s AND -1.4
                    group by lADCP.PositionLat_dec
                    order by lADCP.PositionLat_dec
                    
                    """%(cruise, lat_s), mydb)
                    

            else:
                df_lADCP = sql.read_sql(""" 
            
                
                select lADCP.PositionLat_dec, lADCP.Depth, lADCP.u
                FROM lADCP_gridded_data as lADCP
                WHERE lADCP.Cruise= %s AND lADCP.Depth BETWEEN 0 AND 6000 
                AND lADCP.PositionLat_dec BETWEEN %s AND %s
                
                """%(cruise, lat_s, lat_n), mydb)
                

            
            print("ladcp")
            #print df_lADCP.head(10)
            
            lat = np.array(df_lADCP['positionlat_dec'])
            depth = np.array(df_lADCP['depth']*-1)
            param_ = np.array(df_lADCP['u'])
            
            if measure == 'POC_hFlux':
                xi = np.linspace(min(lat_grid),max(lat_grid),140)
                yi = np.linspace(min(depth_grid),max(depth_grid),160)
                zi_adcp = griddata((lat,depth), param_, (xi[None,:], yi[:,None]), method=grid_method)
            
            
        #### particle data
        
        if item not in ['tricho', "prot", "zoo", 'cops']:
        
            df = sql.read_sql(""" 
            select u.PositionLon_dec, u.PositionLat_dec, u.Depth, u.uvp_%s
            FROM
            (select uvp.PositionLat_dec, uvp.PositionLon_dec, uvp.CTD_filename, uvp.Depth, %s as uvp_%s
            FROM uvp5_bin_data as uvp
            WHERE uvp.Cruise= %s AND uvp.Depth BETWEEN %s AND  %s) as u
            
            WHERE u.PositionLat_dec BETWEEN %s AND %s
            and u.PositionLon_dec BETWEEN %s AND %s
            
            """%(measure_sql, sql_select, measure_sql, cruise, depth_tresh_upper, depth_tresh_lower, lat_s, lat_n, lon_w, lon_e), mydb)
            
            df = df.fillna(0)
        
        else:
        
        
            df = sql.read_sql(""" 
            select max(u.PositionLat_dec) as positionlat_dec, u.Depth - mod(u.Depth::NUMERIC, 20) + 10) as depth, avg(u.uvp_%s) as uvp_%s
            FROM
            (select uvp.PositionLat_dec, uvp.PositionLon_dec, uvp.CTD_filename, uvp.Depth, %s as uvp_%s
            FROM uvp5_bin_data as uvp
            WHERE uvp.Cruise= %s AND uvp.Depth BETWEEN %s AND  %s) as u
            
            WHERE u.PositionLat_dec BETWEEN %s AND %s
            and u.PositionLon_dec BETWEEN %s AND %s
            
            group by u.Depth - mod(u.Depth::NUMERIC, 20) + 10)
            
            """%(measure_sql, measure_sql, sql_select, measure_sql, cruise, depth_tresh_upper, depth_tresh_lower, lat_s, lat_n, lon_w, lon_e), mydb)
            
            df = df.fillna(0)
        
        #####converting the data into numpy arrays
        
        
        if measure not in ["ADCP"]:
            lat = np.array(df['positionlat_dec'])
            lat_distinct = np.unique(lat)
            if measure == 'POC_hFlux':
                param_ = np.array(df['uvp_poc'])
            else:
                param_ = np.array(df[select_dict[measure]])
            depth = dpth(np.array(df['depth']),lat)*-1
        
        
        print(len(param_))
        print(min(depth))
        
        if len(param_) > 0:
            
            print(max(param_))
            
            #### gridding
            
            
            xi = np.linspace(min(lat_grid),max(lat_grid),140)
            yi = np.linspace(min(depth_grid),max(depth_grid),160)
            zi_ = griddata((lat,depth), param_, (xi[None,:], yi[:,None]), method=grid_method)
            
            if measure == 'POC_hFlux':
                zi_ = zi_*zi_adcp
                print("POC_hFlux: " + str(np.nanmean(zi_)))
                
            
            #### subplotting
            
            ax_ = fig.add_axes([plot_position, fig_level, width, height], ylim = (set_ylim_lower,set_ylim_upper), xlim = xlimit_)
            ax_1 = plot2 = plt.contourf(xi,yi,zi_,contour_levels,cmap=cmap_)
            # draw colorbar
            plt.scatter(lat_distinct, np.repeat((set_ylim_upper), len(lat_distinct)), c = "k", marker = "|", s = 20)
            ####plot data points.
            #plt.scatter(lat,depth,marker='.',c='grey',s=0.1,zorder=10)
            
            
            #if cruise_counter == 5:
            #    ax_.set_title(title, fontsize = 6)
                
            ax_.text(-6.5, -50, letter_dict[cruise_counter], color = "k", fontsize = 8, fontweight = "bold")
            
            if cruise_counter not in  [1,4]:
                ax_.axes.get_yaxis().set_visible(False)
            if 'Tara' in c and p_p == 2:
                ax_.axes.get_yaxis().set_visible(True)
            if cruise_counter in [4,5,6]:  
                ax_.axes.get_xaxis().set_visible(False)
            if cruise_counter == 2:
                 plt.xlabel("Latitude")
                
            
            
            #min_depth = min(bottom_topo_depth)
            
            
            if measure == 'ADCP':
                if c in cruises_23W:
                    ax_text(-1, -180, 'EUC', color = "k", fontsize = 8, fontweight = "bold")
                    ax_.text(-1, -500, 'EIC', color = "k", fontsize = 8, fontweight = "bold")
                    ax_.text(6, -500, 'eastward', color = "white", fontsize = 8, fontweight = "bold", rotation = "vertical", zorder = 1)
                    ax_.text(6, -2000, 'westward', color = "white", fontsize = 8, fontweight = "bold", rotation = "vertical")   
                      
            
            
                if c == 'RB_p16n':
                    ax_.fill_between(df_lADCP_min['lat'].tolist(), df_lADCP_min['max_depth'].tolist(), y2 = int(depth_select)*-1, facecolor = 'white', edgecolor="None")
            
            ax_.fill_between(bottom_topo_lat, bottom_topo_depth, y2 = int(depth_select)*-1, facecolor="k", edgecolor="None")
            
            
             
        
    cruise_counter += 1    

cbaxes = fig.add_axes([0.18, 0.065, 0.7, 0.025]) 
cb = plt.colorbar(ax_1, cax = cbaxes, orientation = 'horizontal')
cb.ax.set_xlabel('vertical POC flux (mg C m^-2 d^-1)')



os.chdir("/Users/rkiko/SVN/uvp5_analysis/Results/Manuscripts/Marine_equatorial_snowfall_v2/plots")

fig_name_pdf = "Fig1_5_5_onlyPOCvflux_v1.pdf"

plt.savefig(fig_name_pdf)

#plt.show()

mydb.close()
plt.close()



"""

"""
