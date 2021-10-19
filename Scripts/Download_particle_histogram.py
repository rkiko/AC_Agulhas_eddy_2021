from paruvpy import db_connect
from paruvpy import read_particle_histogram_data_floats
import sys


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    import getpass
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    user_ = getpass.getpass('Username (note, if you want to exit the script Download_particle_histogram, please write "e" as input:')
    if user_=='e':  print('exited from Download_particle_histogram');exit()
    passwd = getpass.getpass('Password (note, if you want to exit the script Download_particle_histogram, please write "e" as input::')
    if passwd == 'e':  print('exited from Download_particle_histogram');exit()
    mydb = db_connect(user_, passwd)
    valid_projectid_list = '(356)'
    df = read_particle_histogram_data_floats(valid_projectid_list, mydb)
    df.to_csv('Ecopart_histogram_data_raw.tsv', sep='\t', index=False)
