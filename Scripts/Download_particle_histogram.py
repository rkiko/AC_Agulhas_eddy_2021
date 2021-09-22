from paruvpy import db_connect
from paruvpy import read_particle_histogram_data
import sys


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    user_ = sys.argv[1]
    passwd = sys.argv[2]

    mydb = db_connect(user_, passwd)
    valid_projectid_list = '(356)'
    df = read_particle_histogram_data(valid_projectid_list, mydb)
    df.to_csv('Ecopart_histogram_data_raw.tsv', sep='\t', index=False)