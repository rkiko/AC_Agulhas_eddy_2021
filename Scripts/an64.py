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
#I load the data created with matlab script
########################################################################################################################

filename1=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_partialFile_an64m.csv" % home).expanduser()
# I read the data
data=pd.read_csv(filename1, sep=',', header=0)

# I define the column names
column_labels = []
for i in range(0,data.shape[1]):
    column_labels.append(data.columns[i])


########################################################################################################################
#I load the nc file to extract the dates
########################################################################################################################

filename2='6903095_Sprof_all.nc'
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
ds = nc.Dataset('%s/%s' % (storedir,filename2))

Date_Num=np.array(ds.variables['JULD'])+matlab_datenum(1950,1,1)
# Date_Num=np.append(Date_Num,[743737])
########################################################################################################################
#I append the values
########################################################################################################################
column_labels.append('Datenum')

data_new = np.zeros((Date_Num.size,data.shape[1]+1 ))
data_new = pd.DataFrame(data=data_new, columns =column_labels)

index_column = column_labels.index('Datenum')
data_new.values[:, index_column] = Date_Num

for i in range(0,data.shape[1]):
    data_new[column_labels[i]] = data[column_labels[i]]

# In case I have more Argo profiles than colocalisation between float and eddy position, I add 0 at the end
data_new.values[data.shape[0]:Date_Num.size,0:data.shape[1]] = 0
# I manually set the position of the 56th profile as outside the eddy because, even if the float is located inside the
# eddy, the temperature and the bbp present an outlier there. Also, before the recalculation of the eddy centroid and
# eddy-float distances by Remi, this profile was located outside the eddy
data_new['sel_insideEddy'][55]=0
########################################################################################################################
#I save
########################################################################################################################
filename_save=Path("%s/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_an64py.csv" % home).expanduser()

data_new.to_csv(filename_save,index=False)

