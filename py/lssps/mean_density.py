import numpy as np
import lssps._lssps as c

def from_grid(grid, xyz, factor):
    """
    Compute nbar at position xyz from given grid
    nbar[i] = factor*grid(x[i]) 

    Args:
      grid (Grid): Grid object of mean number density
      xyz (array): Array of positions; 3 columns for x, y, z
      factor (float): nbar = factor*grid(x)

    Returns:
      nbar (array): Array of nbar at position x
    """

    nbar = np.empty(xyz.shape[0])
    
    c._mean_density_from_grid(grid._grid, xyz, factor, nbar)

    return nbar

class Points:
    def __init__(self, xyz, w=None):
        self._points= c._kdpoints_alloc();
        c._kdpoints_from_array(self._points, xyz, w)

    def __len__(self):
        return c._kdpoints_len(self._points)
        
    @property
    def density(self):
        a = np.zeros(len(self))
        c._kdpoints_density_as_array(self._points, a)
        return a

    @property
    def density_averaged(self):
        return None

class KDTree:
    def __init__(self, points, *, quota=16):
        self._kdtree= c._kdtree_alloc(points._points, quota)

    def estimate_density_adaptive(self, points, knbr):
        if not (knbr > 0):
            raise ValError('knbr must be positive; knbr= %d' % knbr)
        
        c._mean_density_adaptive_estimate(self._kdtree, points._points, knbr)

