from paruvpy import ecopart_connect
from paruvpy import read_particle_histogram_data_floats
import sys


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    import getpass
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    id_project = int(sys.argv[1])
    filename_Download = sys.argv[2]

    user_ = getpass.getpass('Username (note, if you want to exit the script Download_particle_histogram, please write "e" as input:')
    if user_=='e':  print('exited from Download_particle_histogram_%d' % id_project);exit()
    passwd = getpass.getpass('Password (note, if you want to exit the script Download_particle_histogram, please write "e" as input::')
    if passwd == 'e':  print('exited from Download_particle_histogram_%d' % id_project);exit()
    mydb = ecopart_connect(user_, passwd)
    valid_projectid_list = '(%d)' % id_project
    df = read_particle_histogram_data_floats(valid_projectid_list, mydb)
    df.to_csv('%s' % filename_Download,sep='\t', index=False)
