'''
Download raw particle data from ecopart
'''

import pandas as pd
from paruvpy import calculate_mip_and_map_abundance_w_header_from_df_func
from paruvpy import calculate_flux_w_header_from_df_func
from paruvpy import calculate_mip_and_map_poc_cont_w_header_from_df_func
from paruvpy import calculate_respi_diffusion_limited_w_header_from_df_func






if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))
    df = pd.read_csv("Ecopart_processed_data.tsv", sep='\t') #use nrows for testing
    df = calculate_flux_w_header_from_df_func(df)
    df = calculate_mip_and_map_abundance_w_header_from_df_func(df)
    df = calculate_mip_and_map_poc_cont_w_header_from_df_func(df)
    temperature = df["Temperature [degrees Celsius]"]
    oxygen = df["Doxy [micromol/kg]"]

    df = calculate_respi_diffusion_limited_w_header_from_df_func(df, temperature, oxygen, 4)

    df = df[["RAWfilename", "Latitude", "Longitude", "Date_Time", "Pressure [dbar]", "Vol [L] (sampled for this depth bin)", "MiP_abun", "MaP_abun", "Mip_POC_cont_mgC_m3", "Map_POC_cont_mgC_m3", "Flux_mgC_m2", "Respi_nmolO2_l_h"]]

    df.to_csv('Ecopart_mip_map_flux_data.tsv', sep='\t', index=False)

