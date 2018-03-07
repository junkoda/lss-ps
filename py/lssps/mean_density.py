import numpy as np
import lssps._lssps as c

def from_grid(grid, x, factor):
    """
    Compute nbar at position xyz from given grid
    nbar[i] = factor*grid(x[i]) 

    Args:
      grid (Grid): Grid object of mean number density
      x (array): Array of positions; 3 columns xyz
      factor (float): nbar = factor*grid(x)

    Returns:
      nbar (array): Array of nbar at position x
    """

    nbar = np.empty(xyz.shape[0])
    
    c._py_mean_density_from_grid(grid, xyz, factor, nbar)

    return nbar
