import lssps._lssps as c
import numpy as np

import lssps

def _compute_delta_l(grid, indices, n=0, *,
                     grid_moment=None, grid_delta_l=None):
    """
    Compute delta_l (l = 2 or 4)
    common operation for compute_delta2() and compute_delta4()

    Args:
      grid (Grid): input grid of delta_k
      indices (List): list of (moment_index, prefactor) pairs
                      e.g. [(0000, 3.0/2.0),] for (3/2) Q_xxxx

      n (int): compute x^{-n} weighted delta_l for window function
      grid_moment: temporaty grid for moments Q_index
      grid_delta_l (Grid): output grid of delta_l

      grid_moment, grid_delta_l can be None

    Returns:
      grid_delta_l (Grid)
    """

    if grid.mode != 'real-space':
        raise ValueError('grid is not in real space')

    if grid_delta_l is None:
        grid_delta_l = lssps.grid.zeros_like(grid)
    else:
        grid_delta_l.clear()

    if grid_moment is None:
        grid_moment = lssps.grid.zeros_like(grid)

    # use same x0 for grid and grid.shifted
    offset = grid.offset*(grid.boxsize/grid.nc)
    x0 = (grid.x0[0] + offset,
          grid.x0[1] + offset,
          grid.x0[2] + offset)

    for idx in indices:
        # convert str '0000' to list [0, 0, 0, 0]
        idx_array = [int(x) for x in idx[0]]

        # assign Q_index(x)
        c._yamamoto_compute_moment_x(grid._grid, x0, idx_array, n,
                                     grid_moment._grid)

        if grid.shifted is not None:
            c._yamamoto_compute_moment_x(grid.shifted._grid, x0, idx_array, n,
                                         grid_moment.shifted._grid)            

        # FFT to Q_index(k)
        grid_moment.fft()

        # interlacing
        if grid.shifted is not None:
            grid_moment.interlaced = False
            grid_moment.interlace()

        # delta_l += (k_i/k) (k_j/k) ... Q_ij..(k)            
        c._yamamoto_compute_moment_k(grid_moment._grid, idx_array, idx[1],
                                     grid_delta_l._grid)

    return grid_delta_l


def compute_delta0(grid, n=0, *, grid_moment=None, grid1=None):
    """
    Compute delta_1 = [ (k_x/k) Q_x + cyc. ]

    Args:
      grid (Grid): input grid of delta_k
      grid1 (Grid): output grid of delta2_k
      grid_moment (Grid): temporaty grid for moment computation

    grid1, grid_moment can be None

    Returns:
      grid0
    """

    return _compute_delta_l(grid, [], n,
                            grid_moment=grid_moment, grid_delta_l=grid1)

def compute_delta1(grid, n=0, *, grid_moment=None, grid1=None):
    """
    Compute delta_1 = [ (k_x/k) Q_x + cyc. ]

    Args:
      grid (Grid): input grid of delta_k
      grid1 (Grid): output grid of delta2_k
      grid_moment (Grid): temporaty grid for moment computation

    grid1, grid_moment can be None

    Returns:
      grid1
    """

    # Pairs of inddecies (xx, yy, zz, xy, yz, zx) and prefactors
    indices = [('0', 1.0), ('1', 1.0), ('2', 1.0)]

    return _compute_delta_l(grid, indices, n,
                            grid_moment=grid_moment, grid_delta_l=grid1)


def compute_delta2(grid, n=0, *, grid_moment=None, grid2=None):
    """
    Compute delta_2 = (3/2) [ (k_x/k)^2 Q_xx + cyc.
                              + 2 (k_x/k)(k_y/k) Q_xy + cyc. ]

    Args:
      grid (Grid): input grid of delta_k
      grid2 (Grid): output grid of delta2_k
      grid_moment (Grid): temporaty grid for moment computation

    grid2, grid_moment can be None

    Returns:
      grid2
    """

    # Pairs of inddecies (xx, yy, zz, xy, yz, zx) and prefactors
    indices = [('00', 1.5), ('11', 1.5), ('22', 1.5), # Q_xx + cyc.
               ('01', 3.0), ('12', 3.0), ('20', 3.0), # Q_xy + cyc.
              ]
    return _compute_delta_l(grid, indices, n,
                            grid_moment=grid_moment, grid_delta_l=grid2)

def compute_delta3(grid, n=0, *, grid_moment=None, grid3=None):
    """
    Compute delta_3 = 
                      

    Args:
      grid (Grid): input grid of delta_k
      grid3 (Grid): output grid of delta2_k
      grid_moment (Grid): temporaty grid for moment computation

    grid3, grid_moment can be None

    Returns:
      grid3
    """

    # Pairs of inddecies (xxx, yyy, zzz, ...) prefactors
    indices = []

    # (5/2) Q_xxx + cyc.
    fac = 2.5
    indices.extend([('000', fac), ('111', fac), ('222', fac)]) # Q_xxx + cyc.

    # (5/2) 3 Q_xxy + cyc. + Q_xxz + cyc.
    fac = 2.5*3
    indices.extend([('001', fac), ('112', fac), ('220', fac),  # Q_xxy + cyc.
                    ('002', fac), ('110', fac), ('221', fac)]) # Q_xxz + cyc.

    # (5/2) 6 Q_xyz
    fac = 2.5*6
    indices.extend([('123', fac),])

    return _compute_delta_l(grid, indices, n,
                            grid_moment=grid_moment, grid_delta_l=grid3)


def compute_delta4(grid, n=0, *, grid_moment=None, grid4=None):
    """
    Compute delta_4

    Args:
      grid (Grid)       : input grid of delta_k
      grid4 (Grid)      : output grid of delta4_k
      grid_moment (Grid): temporaty grid for moment computation

    Returns:
      grid4
    """

    # Pairs of moment indecies (xxxx, ...) and prefactors
    indices = []

    # (35/8) Q_xxxx + cyc.
    fac = 35.0/8.0
    indices.extend([('0000', fac), ('1111', fac), ('2222', fac)])

    # (35/8) 4 Q_xxxy + cyc.
    fac = 35.0/2.0
    indices.extend([('0001', fac), ('1112', fac), ('2220', fac)])
    indices.extend([('0002', fac), ('1110', fac), ('2221', fac)])

    # (35/8) 6 Q_xxyy + cyc.
    fac = 3.0*35.0/4.0
    indices.extend([('0011', fac), ('1122', fac), ('2200', fac)])

    # (35/8) 12 Q_xxyz + cyc.
    fac = 3.0*35.0/2.0
    indices.extend([('0012', fac), ('1120', fac), ('2201', fac)])

    assert(len(indices) == 15)

    return _compute_delta_l(grid, indices, n, grid_moment=grid_moment,
                            grid_delta_l=grid4)


def compute_yamamoto(grid_delta, kind, *,
                     k_min=0.0, k_max=1.0, dk=0.01,
                     n=0, grid_delta_x=None,
                     subtract_shotnoise=True,
                     correct_mas=True):
    """
    Args:
      grid_delta (Grid)
      kind (str): 's' for Scoccimaro hexadecapole
                  'b' or '4' for Bianchi hexiadecapole
                  '1' for dipole, '3' for tripole
                  monopole and quadrupole are always calculated

      n (int): 1/x^n weight for window function multipoles for wide-angle effect

      grid_delta_x (Grid): grid in real-space when grid_delta is in fourier-space
      grid_moment (Grid): temporary grid for moment grid
      grid2 (Grid): temporary grid for grid2
      grid4 (Grid): temporary grid for grid4

    Note:
      The content of the temporary grids grid_moment, grid2, grid4, are
      arbitrary; they will be overwritten.

    Returns:
      PowerSpectrum object with multipoles
    """

    if grid_delta.mode == 'real-space':
        if grid_delta_x is None:
            grid_delta_x = grid_delta
    else:
        if grid_delta_x is None:
            raise ValueError('real-space grid_delta is not provided')

    assert(grid_delta_x.mode == 'real-space')

    kind = kind.lower()
    if 's' in kind and 'b' in kind:
        raise ValueError('scoccimarro (s) and bianchi (b) cannot specified simultaneously in kind %s' % kind)

    grid2 = compute_delta2(grid_delta_x, n)

    if n == 0:
        grid0 = grid_delta
    else:
        grid0 = compute_delta0(grid_delta_x, n)
        

    if '1' in kind or '3' in kind:
        grid1 = compute_delta1(grid_delta_x, n)

    if '3' in kind:
        grid3 = compute_delta3(grid_delta_x, n)
    
    if 'b' in kind or '4' in kind:
        grid4 = compute_delta4(grid_delta_x, n)

    # Fourier transform delta
    if grid_delta.mode == 'real-space':
        grid_delta.fft()

    if grid_delta.shifted is not None:
        grid_delta.interlace()

    if 'b' in kind or '4' in kind:
        # quadrupole and Bianchi hexaxecapole
        _ps = c._power_spectrum_compute_yamamoto(k_min, k_max, dk,
                                                 grid_delta._grid,
                                                 grid0._grid,
                                                 grid2._grid, grid4._grid,
                                                 subtract_shotnoise,
                                                 correct_mas)        
    else:
        # quadrupole and Scoccimarro hexadecapole
        _ps = c._power_spectrum_compute_yamamoto(k_min, k_max, dk,
                                                 grid_delta._grid,
                                                 grid0._grid,
                                                 grid2._grid, None,
                                                 subtract_shotnoise,
                                                 correct_mas)

    if ('1' in kind) and ('3' in kind):
        # dipole and tripole
        c._power_spectrum_compute_yamamoto_odd(_ps, k_min, k_max, dk,
                                               grid_delta._grid,
                                               grid1._grid,
                                               grid3._grid,
                                               subtract_shotnoise,
                                               correct_mas)
    elif '1' in kind:
        # dipole
        c._power_spectrum_compute_yamamoto_odd(_ps, k_min, k_max, dk,
                                               grid_delta._grid,
                                               grid1._grid, None,
                                               subtract_shotnoise,
                                               correct_mas)

        
    return lssps.PowerSpectrum(_ps)

