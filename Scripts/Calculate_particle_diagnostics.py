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
    import warnings
    import numpy as np
    from pathlib import Path
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    filename_processed = sys.argv[1]
    filename_diagnostics = sys.argv[2]
    id_project = int(filename_processed.split('_')[-1].split('.')[0])

    df = pd.read_csv("%s" % filename_processed, sep='\t')
    ####################################################################################################################
    # Here I calculate the flux
    ####################################################################################################################
    df = calculate_flux_w_header_from_df_func(df)
    df = calculate_flux_w_header_from_df_func(df,eta=0.62,b=66)
    df = calculate_flux_w_header_from_df_func(df,lower_limit=0.02539)
    df = calculate_flux_w_header_from_df_func(df,lower_limit=0.02539,eta=0.62,b=66)
    ####################################################################################################################
    # Here I calculate the MiP and MaP abundances and POC content
    ####################################################################################################################
    df = calculate_mip_and_map_abundance_w_header_from_df_func(df)
    df = calculate_mip_and_map_poc_cont_w_header_from_df_func(df)
    df = calculate_mip_and_map_poc_cont_w_header_from_df_func(df,lower_limit_mip = 0.02539)

    ####################################################################################################################
    # Here I calculate the respiration rates in different ways
    ####################################################################################################################
    temperature = df["Temperature [degrees Celsius]"].copy()
    oxygen = df["Doxy [micromol/kg]"].copy()

    warnings.filterwarnings("ignore", message="A value is trying to be set on a copy of a slice from a DataFrame")
    sel = (temperature == 99999) & (oxygen == 99999)
    temperature[sel] = np.nan
    oxygen[sel] = np.nan

    list_kRemPoc=[0.01300, 0.00300, 0.03100, 0.10000, 0.50000]
    list_lower_limit0=[0.1, 0.02539]
    for lower_limit in list_lower_limit0:
        df = calculate_respi_diffusion_limited_w_header_from_df_func(df, temperature, oxygen, 4,'Kalvelage',lower_limit0 = lower_limit)
        df = calculate_respi_diffusion_limited_w_header_from_df_func(df, temperature, oxygen, 4,'Iversen',lower_limit0 = lower_limit)
        ikRemPoc = list_kRemPoc[0]
        for ikRemPoc in list_kRemPoc:
            df = calculate_respi_diffusion_limited_w_header_from_df_func(df, temperature, oxygen, 4, 'Reminer', kRemPoc=ikRemPoc,lower_limit0 = lower_limit)

    ####################################################################################################################
    # Here I define the list of columns which I save: I take all the columns in which I have the mip, the map, the flux.
    # or the respiration rates
    ####################################################################################################################
    list_columns = ["RAWfilename", "Latitude", "Longitude", "Date_Time", "Pressure [dbar]",
                    "Vol [L] (sampled for this depth bin)"]
                    #, "MiP_abun", "MaP_abun", "Mip_POC_cont_mgC_m3",
                    # "Map_POC_cont_mgC_m3", "Flux_mgC_m2", "Respi_nmolO2_l_h", "Respi_Iversen_nmolO2_l_h",
                    # "Respi_Reminer_kRemPoc0013_nmolO2_l_h", "Respi_Reminer_kRemPoc0003_nmolO2_l_h",
                    # "Respi_Reminer_kRemPoc0031_nmolO2_l_h", "Respi_Reminer_kRemPoc0100_nmolO2_l_h",
                    # "Respi_Reminer_kRemPoc0500_nmolO2_l_h"]
    df_header = list(df.columns)
    for x in df_header:
        if ("Flux_" in x or "Respi_" in x or "Mip_" in x or "Map_" in x or "MiP_" in x or "MaP_" in x) and '_add' not in x:
            list_columns.append(x)


    # For the data coming from a cruise, I add also the cruise and
    # ctd_filename columns. For the uvp6 on equipped on floats, this is not necessary
    list_Cruise_Projects = [231]
    if id_project in list_Cruise_Projects:
        list_columns.insert(1,"Cruise")
        list_columns.insert(1,"CTD_filename")

    df = df[list_columns]
    # I remove the duplicated columns: indeed, the MaP calculation does not change when I extend it to the smaler size
    # classes, because this affects the MiP only
    df = df.loc[:, ~df.columns.duplicated()]

    df.to_csv('%s' % filename_diagnostics, sep='\t', index=False)

