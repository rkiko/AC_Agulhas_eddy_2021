#'CALCULATE OXYGEN CONCENTRATION AT SATURATION
#'Source:  Oxygen solubility in seawater: Better fitting equations.
#'          Garcia & Gordon (1992) L&O 37: 1307-1312.
#'Input:       S = Salinity (pss-78)
#'             T = Temp (deg C)
#'Output:      Oxygen saturation at one atmosphere (umol/kg).

#'DEFINE CONSTANTS, ETC FOR SATURATION CALCULATION

#' The constants used are for units of umol O2/kg.

# O2sat_new(35,20)           = 225.5366
# own_func.O2sat_func(35,20) = 225.53661514454456

def O2sat_func(sal, temp):
    import numpy as np
    A0, A1, A2, A3, A4, A5 = 5.80871, 3.20291, 4.17887, 5.10006, -0.0986643, 3.80369

    B0, B1, B2, B3 = -0.00701577, -0.00770028, -0.0113864, -0.00951519

    C0 = -0.000000275915

    ####### conducting the aou calculation#########
    TS = np.log((298.15 - temp) / (273.15 + temp))
    A = ((((A5 * TS + A4) * TS + A3) * TS + A2) * TS + A1) * TS + A0
    B = ((B3 * TS + B2) * TS + B1) * TS + B0
    O2sat = np.exp(A + sal * (B + sal * C0))
    return O2sat

