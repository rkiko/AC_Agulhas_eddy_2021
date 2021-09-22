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
    test = 0
    if test == 1:
        df = pd.read_csv("Ecopart_histogram_data_raw.tsv", sep='\t', nrows = 200) #use nrows for testing
    else:
        df = pd.read_csv("Ecopart_histogram_data_raw.tsv", sep='\t') #use nrows for testing
    df = process_raw_ecopart_histo_df(df)
    df = calculate_flux_func(df)
    df = calculate_mip_and_map_abundance_func(df)
    df = df[["RAWfilename", "Latitude", "Longitude", "Date_Time", "Pressure [dbar]", "Vol [L] (sampled for this depth bin)", "MiP_abun", "MaP_abun", "Flux_mgC_m2"]]
    if test == 1:
        df.to_csv('Ecopart_mip_map_flux_data_test.tsv', sep='\t', index=False)
    else:
        df.to_csv('Ecopart_mip_map_flux_data.tsv', sep='\t', index=False)

