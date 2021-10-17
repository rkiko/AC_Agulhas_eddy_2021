


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_scripts = Path('~/GIT/AC_Agulhas_eddy_2021/Scripts').expanduser()
    os.chdir(str(path_to_scripts))
    #os.system("python Download_particle_histogram.py")
    os.system("python Convert_raw_particle_data_to_pangaea_style_dataset.py")
    os.system("python Append_BGC_variables.py")
    os.system("python Calculate_particle_diagnostics.py")

