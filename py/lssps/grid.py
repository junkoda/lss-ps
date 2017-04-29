import lssps
import lssps._lssps as c
from numbers import Number
import warnings


class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension

    Methods:
        grid[:]: return grid data as an array
        clear(): reset all data to 0
        fft():   Fast Fourier Transform to Fourier space

    Properties:
        grid.boxsize:
        grid.x0:
        grid.offset:
        grid.shifted: A grid shifted by half grid spacing (for interlacing)
    """

    def __init__(self, nc, boxsize, x0=None, offset=None):
        # _grid is the pointer to a C++ Grid object
        self._grid = c._grid_alloc(nc)
        self.shifted = None
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
        if self.shifted is not None:
            self.shifted.clear()

        return self

    def fft(self):
        """Fast Fourier Transform from real space to Fourier space"""
        c._grid_fft(self._grid)
        if self.shifted is not None:
            self.shifted.fft()

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
        """Length of the cubic box on a side"""
        return c._grid_get_boxsize(self._grid)

    @boxsize.setter
    def boxsize(self, value):
        c._grid_set_boxsize(self._grid, value)
        if self.shifted is not None:
            self.shifted.boxsize = value

    @property
    def x0(self):
        """Coordinate of the box corner"""
        return c._grid_get_x0(self._grid)

    @x0.setter
    def x0(self, x0):
        """Coordinate of the box corner
           grid.x0 = (0.0, 0.0, 0.0)"""
        c._grid_set_x0(self._grid, x0[0], x0[1], x0[2])
        if self.shifted is not None:
            self.shifted.x0 = x0

    @property
    def offset(self):
        """offset of the grid with respect to the box
           in units of grid spacing"""

        return c._grid_get_offset(self._grid)

    @offset.setter
    def offset(self, value):
        c._grid_set_offset(self._grid, value)
        if self.shifted is not None:
            self.shifted.offset = value - 0.5


def zeros(nc, boxsize, x0=None, offset=0.0, *, interlacing=False):
    """Return a new empty grid filled with zeros
    Args:
        nc (int): number of grids per dimension
        boxsize: length of the box on a side
        x0: corner of the box
        offset: offset of the grid points from the box corner
                give a tuple of 2 floats to set both offsets for
                the main and shifted grid
        interlacing: attach shifted grid which is shifted by half gridspacing
    """

    grid = Grid(nc, boxsize, x0, offset)
    if interlacing:
        grid.shifted = Grid(nc, x0, boxsize, offset - 0.5)

    grid.clear()

    return grid


def compute_fluctuation(grid_data, grid_rand=None):
    """Compute density fluctuation data = data - rand

    Args:
        data: a tuple of data grids
        rand: a tuple of random grids (may be None)

    Returns:
        fluctutation grid (= grid_data)
        data grids are modified and become fluctuation grid
    """

    if grid_rand is None:
        c._grid_compute_fluctuation_homogeneous(grid_data._grid)
    else:
        c._grid_compute_fluctuation(grid_data._grid, grid_rand._grid)

    if grid_data.shifted is not None:
        if grid_rand is None:
            compute_fluctuation(grid_data.shifted)
        else:
            compute_fluctuation(grid_data.shifted, grid_rand.shifted)

    return grid_data
