import lssps._lssps as c


class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension"""

    def __init__(self, nc):
        # _grid is the C++ Grid object
        self._grid = c._grid_alloc(nc)

    def __getitem__(self, index):
        """grid[ix, iy, iz]: value at grid point (ix, iy, iz)
           grid[:]         : whole grid as a numpy.array"""
        if self.mode == 'real-space':
            return c._grid_fx_asarray(self._grid)[index]
        elif self.mode == 'fourier-space':
            return c._grid_fk_asarray(self._grid)[index]

    def clear(self):
        """Reset the grid with zeros"""
        c._grid_clear(self._grid)
        return self

    def fft(self):
        """Fast Fourier Transform from real space to Fourier space"""
        c._grid_fft(self._grid)
        return self

    @property
    def mode(self):
        """Current mode of the grid:
        'real-space' or 'fourier-space'
        """
        return c._grid_mode(self._grid)

    @property
    def nc(self):
        """Number of grid points per dimension"""
        return c._grid_nc(self._grid)


def zeros(nc):
    """Return a new empty grid filled with zeros"""
    grid = Grid(nc)
    grid.clear()

    return grid
