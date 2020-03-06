import numpy as np


def is_bgs_goodz(z, zerr, ztruth):
    ##  The random error on BGS redshifts in a Gaussian core shall be less than
    ##  sigma_z = 0.0005 * (1. + z), or 150 km/s.
    is_good = zerr < 0.0005 * (1. + ztruth)  

    ##  Systematic inaccuracy in the mean BG redshift shall be less than
    ##  $\Delta\bar{z} = 0.0004(1. + z)$ (equivalent to 120 km/s).
    is_good = is_good & 
    
    return  is_good
