import lssps
import lssps._lssps as c
from numbers import Number
import warnings
import h5py

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
        self.interlacing = None

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
        self.mode = 'real-space'
        
        if self.shifted is not None:
            self.shifted.clear()

        return self

    def compute_fluctuation(self, grid_rand=None):
        """Compute density fluctuation data -= rand

        Args:
           grid_rand: the grid of randoms  (may be None)
        """

        if grid_rand is None:
            c._grid_compute_fluctuation_homogeneous(self._grid)
        else:
            c._grid_compute_fluctuation(self._grid, grid_rand._grid)

        if self.shifted is not None:
            if grid_rand is None:
                self.shifted.compute_fluctuation()
            else:
                self.shifted.compute_fluctuation(grid_rand.shifted)

        return self

    def interlace(self):
        """perform interlacing"""

        if self.shifted is None:
            raise TypeError('grid.shifted is required for interlacing')
        if self.mode != 'fourier-space':
            raise RuntimeError('grid must be in Fourier space for interlacing')

        c._interlacing(self._grid, self.shifted._grid)
        self.interlacing = 'done'

        return self

    def fft(self):
        """Fast Fourier Transform from real space to Fourier space"""
        c._grid_fft(self._grid)
        if self.shifted is not None:
            self.shifted.fft()

        return self

    def save_h5(self, filename):
        """Save grid to an HDF5 file
        
        """
        f = h5py.File(filename, 'w')
        
        f['nc'] = self.nc
        f['boxsize'] = self.boxsize
        f['x0_box'] = self.x0
        f['offset'] = self.offset
        sums = self.sums

        f['sums'] = (sums[0], sums[1], sums[2])
        f['np'] = sums[3]
        f['n_mas'] = self.n_mas
        
        if self.mode == 'real-space':
            f['fx'] = self[:]
        elif self.mode == 'fourier-space':
            f['fk'] = self[:]

        if self.shifted is not None:
            g = f.create_group('shifted')
            g['offset'] = self.shifted.offset
            sums = self.shifted.sums
            g['sums'] = (sums[0], sums[1], sums[2])
            g['np'] = sums[3]
            
            if self.mode == 'real-space':
                g['fx'] = self.shifted[:]
            elif self.mode == 'fourier-space':
                g['fk'] = self.shifted[:]

        f.close()
    # def assign_density(cat, mas):
    #     """Assign density from particles in file or array"""

    #     n_mas = lssps._mass_assignment_scheme[mas.upper()]
        
    #     # From CatalogueFile
    #     if isinstance(catalogue_file, CatalogueFile):
    #         c._mass_assignment(cat._f, x0, boxsize,
    #                            n_mas, grid._grid, None)


    #     if interlacing:
    #         grid_shifted = lssps.grid.zeros(nc)
    #         c._mass_assignment(catalogue_file._f, x0, boxsize,
    #                            n_mas, grid._grid, grid_shifted._grid)

    #     return (grid, grid_shifted)

    #if isinstance(cat, str):\
        
    # catalogue is the filename
    #    cat = loadtxt(cat)





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

    @property
    def sums(self):
        """Length of the cubic box on a side"""
        ss = c._grid_get_sums(self._grid)
        return c._grid_get_sums(self._grid)

    @sums.setter
    def sums(self, values):
        c._grid_set_sums(self._grid, values[0], values[1], values[2], values[3])

    @property
    def n_mas(self):
        """Length of the cubic box on a side"""
        return c._grid_get_nmas(self._grid)

    @n_mas.setter
    def n_mas(self, value):
        c._grid_set_nmas(self._grid, value)




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
        grid.shifted = Grid(nc, boxsize, x0, offset - 0.5)

    grid.clear()

    return grid



def load_h5(filename):
    """
    Args:
        filename (str): file name of an HDF file
    """
    
    if isinstance(filename, str):
        f = h5py.File(filename, 'r')
    else:
        f = filename
        
    nc = f['nc'][()]
    boxsize = f['boxsize'][()]
    x0 = f['x0_box'][:]
    offset = f['offset'][()]

    grid = Grid(nc, boxsize, x0, offset)

    sums = f['sums'][:]
    n = int(f['np'][()])
    n_mas = f['n_mas'][()]

    grid.sums = (sums[0], sums[1], sums[2], n)
    grid.n_mas = n_mas

    a = f['fx'][:]
    if 'fx' in f:
        c._grid_load_fx_from_array(grid._grid, a)
    else:
        raise RuntimeError('fx does not exists in HDF5 file')

    if 'shifted' in f:
        g = f['shifted']
        offset_shifted = g['offset'][()]
        grid.shifted = Grid(nc, boxsize, x0, offset_shifted)

        sums = g['sums'][:]
        n = int(g['np'][()])

        grid.shifted.sums = (sums[0], sums[1], sums[2], n)
        grid.shifted.n_mas = n_mas

        a_shifted = g['fx'][:]
        c._grid_load_fx_from_array(grid.shifted._grid, a_shifted)
    
    f.close()

    return grid


