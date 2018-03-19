import lssps._lssps as c

def gaussian_smoothing(grid, sigma):
    """
    Apply Gaussian smoothing (
    Args:
      grid (Grid): grid of white noise
      sigma (float): smoothing length
    """

    if grid.mode != 'fourier-space':
        raise ValueError('grid must be in fourier-space')

    c._etc_gaussian_smoothing(grid._grid, sigma)
    
