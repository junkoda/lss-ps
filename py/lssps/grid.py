import lssps
import lssps._lssps as c


class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension"""

    def __init__(self, nc, boxsize, x0=None, *, offset=None):
        # _grid is the pointer to a C++ Grid object
        self._grid = c._grid_alloc(nc)
        self.boxsize = boxsize

        if x0 is not None:
            self.x0 = x0

        if offset is not None:
            self.offset = offset

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

    @property
    def boxsize(self):
        """length of the cubic box on a side"""
        return c._grid_get_boxsize(self._grid)

    @boxsize.setter
    def boxsize(self, value):
        c._grid_set_boxsize(self._grid, value)

    @property
    def x0(self):
        """corner coordinate of the box"""
        return c._grid_get_x0(self._grid)

    @x0.setter
    def x0(self, x0):
        """corner coordinate of the box"""
        return c._grid_set_x0(self._grid, x0[0], x0[1], x0[2])

    @property
    def offset(self):
        """offset of the grid with respect to the box in units of grid spacing"""
        return c._grid_get_offset(self._grid)

    @offset.setter
    def offset(self, value):
        return c._grid_set_offset(self._grid, value)


def zeros(nc):
    """Return a new empty grid filled with zeros"""
    grid = Grid(nc)
    grid.clear()

    return grid

def compute_fluctuation(grid_data, grid_rand = None):
    """Compute density fluctuation data = data - rand

    Args:
        data: a tuple of data grids
        rand: a tuple of random grids (may be None)

    Returns:
        a tuple of fluctutation grid

    Note:
        data grids are modified and become fluctuation grids
    """

    
    if isinstance(grid_data, lssps.Grid):
        grid_data = (grid_data,)

    if isinstance(grid_rand, lssps.Grid):
        grid_rand = (grid_rand,)
    
    n = len(grid_data)
    assert(n == 1 or n == 2)
    
    if grid_rand is None:
        for i in range(n):
            c._grid_compute_fluctuation_homogeneous(grid_data[i]._grid)
    else:
        assert(len(grid_rand) == n)
        for i in range(n):
            c._grid_compute_fluctuation(grid_data[i]._grid, grid_rand[i]._grid)

    if len(grid_data) == 1:
        return grid_data[0]
    
    return grid_data
