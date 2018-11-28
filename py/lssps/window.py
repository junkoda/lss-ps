"""
Window function module
"""

import numpy as np
import lssps._lssps as c

def compute_grid3d(grid, *, pk_fac=None, shot_noise=None):
    """
    Compute 3D window function grid |W(k)|^2
    """

    if grid.mode == 'real-space':
        grid.fft()

    if grid.shifted and not grid.interlaced:
        grid.interlace()

    grid.correct_mas()

    #pk_fac = #nbar_mid*grid_rand.w2_sum/alpha
    #shot_noise = pk_fac*grid_rand.w2_sum

    if pk_fac is None:
        pk_fac = grid.pk_normalisation
    if shot_noise is None:
        shot_noise = grid.shot_noise
    
    g = grid[:]
    w2_grid = pk_fac*(g*np.conj(g)) - shot_noise

    return w2_grid

def bessel_transform(r, k, f, l):
    """
    Compute spherical-Bessel transform
      1/(2 pi^2) \int k^2 dk j_l(kr) f(k)

    Args:
      r (array): output r
      k (array): k of f(k)
      f (array): function f
      l (int):   degree of Bessel function j_l
      
    """
    result = np.empty_like(r)
    
    c._window_cquag_bessel_transform(r, k, f, l, result)
    
    return result

