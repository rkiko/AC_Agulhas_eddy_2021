import numpy as np
import datetime
from datetime import date

def matlab_datenum(*argv):
    if argv.__len__()==1:
        x=argv[0]
        if (type(x) is list)|(type(x) is tuple):
            x = np.array(x)

        if x.size==3:
            year,month,day,hour,minute,seconds=x[0],x[1],x[2],0,0,0
        elif x.size==6:
            year,month,day,hour,minute,seconds=x[0],x[1],x[2],x[3],x[4],x[5]
        else:
            raise ValueError("The input array or list should have either 3 or 6 elements")

    elif argv.__len__()==3:
        year,month,day,hour,minute,seconds=argv[0],argv[1],argv[2],0,0,0
    elif argv.__len__()==6:
        year,month,day,hour,minute,seconds=argv[0],argv[1],argv[2],argv[3],argv[4],argv[5]
    else:
        raise ValueError("The input arguments should be either 3 or 6 elements, either a 3 or 6 element list/array")

    date_reference = datetime.datetime.strptime("%02d/%02d/%d" % (day,month,year), "%d/%m/%Y")
    date_reference_datenum = date.toordinal(date_reference) + 366 +hour/24+minute/24/60+seconds/86400
    return date_reference_datenum

if __name__ == "__main__":
    year=1950
    month=1
    day=1
    hour=12
    x=[year, month, day, hour, 0, 0]
    result = matlab_datenum(x)
    print(result)
    result = matlab_datenum(year, month, day, hour, 0, 0)
    print(result)
    result = matlab_datenum(year, month, day, 0, 0, 0)
    print(result)

