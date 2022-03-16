import numpy as np
import datetime
from datetime import date

def matlab_datevec(x):

    fracx=x-np.floor(x)
    y = date.fromordinal(int(np.floor(x)-366))
    hour=np.floor(fracx*24)
    hourfrac=fracx*24-hour
    minute=np.floor(hourfrac*60)
    minutefrac=hourfrac*60-minute
    seconds=minutefrac*60
    result=np.array([y.year,y.month,y.day,hour,minute,seconds])
    return result

if __name__ == "__main__":
    x=734000.5567
    result = matlab_datevec(x)
    print(result)

