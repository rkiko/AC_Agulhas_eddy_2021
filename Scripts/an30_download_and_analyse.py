


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    path_to_scripts = Path('~/GIT/AC_Agulhas_eddy_2021/Scripts').expanduser()
    os.chdir(str(path_to_scripts))
    home = str(Path.home())
    sys.path.insert(0, "%s/GIT/Lagrangian_uvp/Scripts" % home)
    # from Download_Ecopart import Download_Ecopart
    id_project=231
    filename0='Ecopart'

    filename_Download='%s_histogram_data_raw_%d.tsv' % (filename0,id_project)
    os.system("python Download_particle_histogram.py %d %s" % (id_project,filename_Download))
    print('Download of project id %d finished, starting conversion of raw data' % id_project)

    filename_processed = '%s_processed_data_%d.tsv' % (filename0,id_project)
    os.system("python Convert_raw_particle_data_to_pangaea_style_dataset.py %s %s" % (filename_Download,filename_processed))
    print('Conversion of raw data of project id %d finished ' % id_project)

    os.system("python Append_CTD_data.py %s" % (filename_processed))
    print('Biogeochemical variables attached to processed data of project id %d ' % id_project)

    filename_diagnostics='%s_diagnostics_data_%d.tsc' % (filename0,id_project)
    os.system("python Calculate_particle_diagnostics.py %s %s" % (filename_processed,filename_diagnostics))
    print('Calculation of diagnostics of project id %d finished' % id_project)

    os.system("python Append_CTD_data.py %s" % (filename_diagnostics))
    print('Biogeochemical variables attached to diagnostics of project id %d ' % id_project)
