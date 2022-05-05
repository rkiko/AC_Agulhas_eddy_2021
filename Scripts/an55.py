import numpy as np
import pandas as pd
import os
from pathlib import Path
home = str(Path.home())
storedir='%s/GIT/AC_Agulhas_eddy_2021/Data' % home
sys.path.insert(0, "%s/GIT/AC_Agulhas_eddy_2021/Scripts" % home)
path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
os.chdir(str(path_to_data))

#######################################################################
# Parameters
#######################################################################
id_project = 356
filename0 = 'Ecopart'
filename_processed = '%s_processed_data_%d.tsv' % (filename0,id_project)
filename_diagnostics = '%s_diagnostics_data_%d.tsv' % (filename0, id_project)

df = pd.read_csv("%s" % filename_processed, sep='\t')
from paruvpy import calculate_mip_and_map_flux_w_header_from_df_func

df = calculate_mip_and_map_flux_w_header_from_df_func(df)
df = calculate_mip_and_map_flux_w_header_from_df_func(df, eta=0.62, b=66)
df = calculate_mip_and_map_flux_w_header_from_df_func(df, lower_limit_mip=0.02539)
df = calculate_mip_and_map_flux_w_header_from_df_func(df, lower_limit_mip=0.02539, eta=0.62, b=66)

# I remove the duplicated columns: indeed, the MaP calculation does not change when I extend it to the smaler size
# classes, because this affects the MiP only
df = df.loc[:, ~df.columns.duplicated()]

