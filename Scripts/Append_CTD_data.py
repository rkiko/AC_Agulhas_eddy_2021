from paruvpy import include_CTD_data


if __name__ == '__main__':
    import os, sys
    from pathlib import Path
    import pandas as pd
    path_to_data = Path('~/GIT/AC_Agulhas_eddy_2021/Data').expanduser()
    os.chdir(str(path_to_data))

    filename = sys.argv[1]
    data = pd.read_csv(filename, sep='\t', header=0)
    data = include_CTD_data(data)
    data.to_csv('%s' % filename, sep='\t', index=False)

