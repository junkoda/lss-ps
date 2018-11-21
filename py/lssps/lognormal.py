"""
Generate Gaussian and lognormal grid
"""

import lssps._lssps as c

def convert_P_to_delta_k(grid, seed, *, fix_amplitude=False):
    """
    Generate random Gaussian delta(k) from P(k) grid

    seed (int): random seed
    fix_amplitude (bool): fix amplitude, only random phase (not Gaussian)
    """

    if not isinstance(seed, int):
        raise TypeError('seed must be an integer: {}'.format(type(seed)))
    
    c._lognormal_convert_P_to_delta_k(grid, seed, int(fix_amplitude))

    return grid
