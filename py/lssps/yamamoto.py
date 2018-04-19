import lssps._lssps as c
import numpy as np

import lssps

def _compute_delta_l(grid, indices, *, grid_moment=None, grid_delta_l=None):
    """
    Compute delta_l (l = 2 or 4)
    common operation for compute_delta2() and compute_delta4()

    Args:
      grid (Grid): input grid of delta_k
      indices (List): list of (moment_index, prefactor) pairs
                      e.g. [(0000, 3.0/2.0),] for (3/2) Q_xxxx
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
        c._yamamoto_compute_moment_x(grid._grid, x0, idx_array,
                                     grid_moment._grid)

        if grid.shifted is not None:
            c._yamamoto_compute_moment_x(grid.shifted._grid, x0,
                                         idx_array,
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

    
def compute_delta2(grid, *, grid_moment=None, grid2=None):
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
    #assert(len(indices) == 6)

    return _compute_delta_l(grid, indices,
                            grid_moment=grid_moment, grid_delta_l=grid2)


def compute_delta4(grid, *, grid_moment=None, grid4=None):
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

    return _compute_delta_l(grid, indices, grid_moment=grid_moment,
                            grid_delta_l=grid4)


def compute_yamamoto(grid_delta, kind, *,
                     k_min=0.0, k_max=1.0, dk=0.01, nmu=0,
                     subtract_shotnoise=True,
                     correct_mas=True,
                     grid_moment=None, grid2=None, grid4=None):
    """
    Args:
      grid_delta (Grid)
      kind (str): 'scoccimarro' or 'bianchi'

      grid_moment (Grid): temporary grid for moment grid
      grid2 (Grid): temporary grid for grid2
      grid4 (Grid): temporary grid for grid4

    Note:
      The content of the temporary grids grid_moment, grid2, grid4, are
      arbitrary; they will be overwritten.

    Returns:
      PowerSpectrum object with multipoles
    """

    if grid_delta.mode != 'real-space':
        raise ValueError('grid is not in real space')

    grid2 = compute_delta2(grid_delta, grid_moment=grid_moment, grid2=grid2)

    if kind.lower() == 'bianchi':
        grid4 = compute_delta4(grid_delta, grid_moment=grid_moment,
                               grid4=grid4)

    # Fourier transform delta
    grid_delta.fft()

    if grid_delta.shifted is not None:
        grid_delta.interlace()


    if kind.lower() == 'scoccimarro':
        _ps = c._power_spectrum_compute_yamamoto(k_min, k_max, dk,
                                                 grid_delta._grid,
                                                 grid2._grid, None,
                                                 subtract_shotnoise,
                                                 correct_mas)
    elif kind.lower() == 'bianchi':
        _ps = c._power_spectrum_compute_yamamoto(k_min, k_max, dk,
                                                 grid_delta._grid,
                                                 grid2._grid, grid4._grid,
                                                 subtract_shotnoise,
                                                 correct_mas)
    else:
        raise ValError('Unknown Yamamoto estimator; must be scoccimarro or bianchi: %s' % kind)

    return lssps.PowerSpectrum(_ps)

