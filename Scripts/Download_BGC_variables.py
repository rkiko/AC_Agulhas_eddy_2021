'''
Download the BGC variables from the Coriolis data center via ftp
'''
from paruvpy import download_Coriolis_data


if __name__ == '__main__':
    import os
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))
    inp=input('To continue, type any key, otherwise, type "e" to exit:')
    if inp=='e': print('exited from Download_BGC_variables');exit()
    valid_BGCArgo_id_list=['6903095']
    download_Coriolis_data(valid_BGCArgo_id_list)
