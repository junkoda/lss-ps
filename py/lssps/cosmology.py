"""
Module for cosmological functions.

D(omega_m, *, a=None, z=None): linear growth factor
f(omega_m, *, a=None, z=None): linear growth rate
H_factor(omega_m, *, a=None, z=None): Hubble parameter / 100 h

compute_comoving_distance(z): compute comoving distance from redshift
compute_redshift(d): compute redshift from comoving distance

"""

import lssps._lssps as c
import numpy as np
import math

def init(omega_m, z_max, n=1001):
    c._cosmology_init(omega_m, z_max, n)

def compute_redshift(d, *, d_min=None, d_max=None, n_iterp=1001):
    """
    Args:
       d (array): array of comoving distances

    Returns:
       z (array): redshifts that correspond to d
    """
    
    z = np.empty_like(d)
    c._cosmology_distance_redshift(d, z)

    return z

def compute_comoving_distance(z, *, z_min=None, z_max=None, n_iterp=1001):
    """
    TODO: change interface, do _init automatically

    Args:
       z (float): redshift
    Returns:
       d (float): comoving distance
    """

    return c._cosmology_compute_comoving_distance(z)

def _get_a(a, z):
    # return scale factor a from either a or redshift z
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

    Provide either a or z.
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

    Provide either a or z.
    """

    a = _get_a(a, z)

    return c._cosmology_growth_f(a, omega_m)

def H_factor(omega_m, *, a=None, z=None):
    """
    Compute redshift-dependent factor in Hubble parameter,
      H(a) = 100 h H_factor [km/s/Mpc],
    which is,
      sqrt(omega_m*a**-3 + (1.0 - omega_m)).

    Args:
       omega_m (float): Omega_m0
       a (float): scale factor
       z (float): redshift

    Provide either a or z.
    """

    a = _get_a(a, z)
    omega_l = 1.0 - omega_m
    
    return math.sqrt(omega_m*a**-3 + omega_l)
