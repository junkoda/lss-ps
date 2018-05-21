import numpy as np
import lssps._lssps as c

def from_grid(grid, xyz, factor, *, kind='NGP'):
    """
    Compute nbar at position xyz from given grid
    nbar[i] = factor*grid(x[i]) 

    Args:
      grid (Grid): Grid object of mean number density
      xyz (array): Array of positions; 3 columns for x, y, z
      factor (float): nbar = factor*grid(x)
    
      kind='TSC' (str): interpolation scheme 'NGP' or 'TSC'

    Returns:
      nbar (array): Array of nbar at position x
    """

    nbar = np.empty(xyz.shape[0])

    if kind.upper() == 'NGP':
        n_interp = 0
    elif kind.upper() == 'CIC':
        n_interp = 1
    elif kind.upper() == 'TSC':
        n_interp = 2
    else:
        raise ValueError('Unknown `kind` for mean_density.from_grid: %s. Must be NGP, CIC or TSC.' % kind)
    
    c._mean_density_from_grid(grid._grid, xyz, factor, n_interp, nbar)

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
    def __init__(self, points, knbr, *, quota=16):
        if not (knbr > 0):
            raise ValError('knbr must be positive; knbr= %d' % knbr)

        self._kdtree= c._kdtree_alloc(points._points, quota)
        c._kdtree_compute_rk(self._kdtree, knbr)
        c._kdtree_update_node_statistics(self._kdtree)

    def estimate_density_adaptive(self, xyz):
        """
        xyz (array [n, 3]): positions for density estiamtion
        """
        a = np.zeros(xyz.shape[0])
        c._mean_density_adaptive_estimate(self._kdtree, xyz, a)

        return a

    def estimate_density_averaged(self, xyz, knbr):
        """
        Compute the average of approximately k neighbors around points xyz

        Args:
          xyz (array): array of positions shape (n_particles, 3)
          knbar (int): number of particles avraging over
        """
        if not (knbr > 0):
            raise ValError('knbr must be positive; knbr= %d' % knbr)
        
        a = np.zeros(xyz.shape[0])
        c._mean_density_average_estimate(self._kdtree, xyz, knbr, a)

        return a
