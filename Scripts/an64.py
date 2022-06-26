import numpy as np
import pandas as pd
import os,sys
import netCDF4 as nc
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory
from matlab_datevec import matlab_datevec
from matlab_datenum import matlab_datenum



########################################################################################################################
#I load the matlab data
########################################################################################################################

filename1=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius.csv" % home).expanduser()
# I read the data
data=pd.read_csv(filename1, sep=',', header=0)

# I define the column names
column_labels = []
for i in range(0,data.shape[1]):
    column_labels.append(data.columns[i])

column_labels.append('Datenum')

data_new = np.zeros((data.shape[0],data.shape[1]+1 ))
data_new = pd.DataFrame(data=data_new, columns =column_labels)

########################################################################################################################
#I load the nc file to extract the dates
########################################################################################################################

filename2='6903095_Sprof_all.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename2))

Date_Num=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
Date_Num=Date_Num[0:data.shape[0]]

########################################################################################################################
#I append the values
########################################################################################################################

index_column = column_labels.index('Datenum')
data_new.values[:, index_column] = Date_Num

for i in range(0,data.shape[1]):
    data_new[column_labels[i]] = data[column_labels[i]]

########################################################################################################################
#I save
########################################################################################################################
filename_save=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/data_an64py.csv" % home).expanduser()

data_new.to_csv(filename_save,index=False)

