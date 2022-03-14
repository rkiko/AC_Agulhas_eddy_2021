import os
from pathlib import Path
home = str(Path.home())
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
os.chdir('%s/GIT/AC_Agulhas_eddy_2021/Scripts/' % home) #changes directory

year=2021
day0=103 #2021,4,13
dayf=266 #2021,9,23

download_webiste='https://oceandata.sci.gsfc.nasa.gov/ob/getfile'
filename0='A'
filename2 = 'L3m_DAY_ZLEE_Zeu_lee_9km.nc'
appkey='a69b30815ad5f6246b09769c97b6cd244f751f6d'

index_day=day0
for index_day in range(day0,dayf+1):
    filename = '%s%d%03d.%s' % (filename0,year,index_day,filename2)
    url_filename_download = '%s/%s?appkey=%s' % (download_webiste, filename,appkey)
    os.system("python obdaac_download.py -v %s" % (url_filename_download))
    os.system("mv %s ../Data/an35/" % (filename))



