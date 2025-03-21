'''
Download raw particle data from ecopart
'''

import pandas as pd
from paruvpy import process_raw_ecopart_histo_df, calculate_mip_and_map_abundance_func
from paruvpy import calculate_sheldon_slope_func, calculate_nbss_slope_func, calculate_flux_func
from paruvpy import calculate_mip_and_map_biovolume_func






if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    filename_input = sys.argv[1]
    filename_output = sys.argv[2]

    test = 0
    if test == 1:
        df = pd.read_csv("%s" % filename_input, sep='\t', nrows=200)  # use nrows for testing
    else:
        df = pd.read_csv("%s" % filename_input, sep='\t') #use nrows for testing
    df = process_raw_ecopart_histo_df(df)
    df.to_csv('%s' % filename_output, sep='\t', index=False)


