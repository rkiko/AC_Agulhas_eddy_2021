'''
Include the BGC variables downloaded from the Coriolis data center via ftp
'''

import pandas as pd
from arguvpy import include_Coriolis_data


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    filename_diagnostics = sys.argv[1]
    id_BGCArgo = str(sys.argv[2])
    df = pd.read_csv(filename_diagnostics, sep='\t')
    df = include_Coriolis_data(df,id_BGCArgo)
    df.to_csv(filename_diagnostics, sep='\t', index=False)


