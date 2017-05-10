import lssps._lssps as c
import numpy as np

def init(omega_m, z_max, n=1001):
    c._cosmology_init(omega_m, z_max, n)

def redshift(d):
    """
    Args:
       d (array): array of comoving distance

    Returns:
       z (array): redshifts that correspond to d
    """
    
    z = np.empty_like(d)
    c._cosmology_distance_redshift(d, z)

    return z
