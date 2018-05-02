import lssps._lssps as c
import numpy as np
import math

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

def compute_comoving_distance(z):
    """
    Args:
       z (float): redshift
    Returns:
       d (float): comoving distance
    """

    return c._cosmology_compute_comoving_distance(z)

def _get_a(a, z):
    if a is None and z is None:
        raise ValueError('Need to give either a or z')

    if z is not None:
        if a is not None:
            raise ValueError('Cannot set both a and z')
        a = 1.0/(1.0 + z)

    return a

def D(omega_m, *, a=None, z=None):
    """
    Compute linear growth factor D(z)

    Args:
       omega_m (float): Omega_m0
       a (float): scale factor
       z (float): redshift
    """

    a = _get_a(a, z)

    return c._cosmology_growth_D(a, omega_m)


def f(omega_m, *, a=None, z=None):
    """
    Compute linear growth rate f(z)

    Args:
       omega_m (float): Omega_m0
       a (float): scale factor
       z (float): redshift    
    """

    a = _get_a(a, z)

    return c._cosmology_growth_f(a, omega_m)

def H_factor(omega_m, *, a=None, z=None):
    """
    Compute redshift-dependent factor in Hubble parameter
    sqrt(omega_m*a**-3 + (1.0 - omega_m))
    H(a) = 100 h H_factor [km/s/Mpc]
    """

    a = _get_a(a, z)
    omega_l = 1.0 - omega_m
    
    return math.sqrt(omega_m*a**-3 + omega_l)
