'''
Download raw particle data from ecopart
'''

from paruvpy import extrapol_small_sizes
import pandas as pd






if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    filename_processed = sys.argv[1]

    df = pd.read_csv(filename_processed, sep='\t')
    df = extrapol_small_sizes(df)
    df.to_csv(filename_processed, sep='\t', index=False)



