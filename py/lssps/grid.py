"""
Grid: class for discritised field in real or Fourier space.

zeros(nc, boxsize, x0=None, offset=0.0, *, interlacing=True)
zeros_like(grid)

empty(nc, boxsize, x0=None, offset=0.0, *, interlacing=True)
empty_like(grid)
"""

import lssps
import lssps._lssps as c
from numbers import Number
import warnings

class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension
    Grid(nc, boxsize, *, x0, offset)

    Methods:
        grid[:]: return grid data as an array
        clear(): reset all data to 0
        fft():   Fast Fourier Transform (FFT) to Fourier space
        fft_inverse(): Inverse FFT from Fourier space to real space
        assign_density(xyz, *, weight=None, nbar=None, mas='CIC')
        compute_fluctuation(grid_rand=None)

    Properties:
        grid.boxsize:
        grid.x0:
        grid.offset: Position of grid points in cell in units of grid spacing
                     (0 <= offset < 1)
        grid.shifted: A grid shifted by half grid spacing for interlacing
        grid.interlaced: True if interlace() has been called. 
                         Reset to None by clear()
    """

    def __init__(self, nc, boxsize, x0=None, offset=None):
        # _grid is the pointer to a C++ Grid object wrapped to _Grid for Python
        self._grid = c._grid_alloc(nc)
        self.shifted = None
        self.boxsize = boxsize
        self.interlaced = False

        if x0 is not None:
            self.x0 = x0   # activates x0.setter

        if offset is not None:
            self.offset = offset

    def __getitem__(self, index):
        """grid[ix, iy, iz]: value at grid point (ix, iy, iz)
           grid[:]         : whole grid as a numpy.array"""
        if self.mode == 'real-space':
            return c._grid_fx_asarray(self._grid)[index]
        elif self.mode == 'fourier-space':
            return c._grid_fk_asarray(self._grid)[index]

    def __setitem__(self, index, value):
        """Set the grid value grid[ix, iy, iz]"""
        if self.mode == 'real-space':
            c._grid_fx_asarray(self._grid)[index] = value
        elif self.mode == 'fourier-space':
            c._grid_fk_asarray(self._grid)[index] = value

    def clear(self):
        """Reset the grid with zeros"""
        c._grid_clear(self._grid)
        self.mode = 'real-space'
        self.interlaced = False
        self.mas_corrected = False
        
        if self.shifted is not None:
            self.shifted.clear()

        return self

    def copy(self):
        """
        Create a copy of this grid.

        Returns:
          a new grid with the same content
        """
        grid = Grid(self.nc, self.boxsize)
        c._grid_copy(self._grid, grid._grid)
        return grid

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
            raise RuntimeError('grid.shifted is required for interlacing')
        if self.mode != 'fourier-space':
            raise RuntimeError('grid must be in Fourier space for interlacing')

        if self.interlaced == False:
            c._interlacing(self._grid, self.shifted._grid)
            self.interlaced = True

        return self

    def fft(self, *, normalise=False):
        """
        Fast Fourier Transform from real space to Fourier space
        
          sum f(x) exp(-ikx)
        
        When normalise=True,
          sum f(x) exp(-ikx) dV
        which is similar to Fourier integral
        """
        if self.mode != 'real-space':
            raise RuntimeError('grid is not in real space: %s' % self.mode)

        c._grid_fft(self._grid)
        if self.shifted is not None:
            self.shifted.fft(normalise=normalise)

        if normalise:
            dV = (self.boxsize/self.nc)**3
            self[:] *= dV

        return self

    def fft_inverse(self, *, normalise=False):
        """
        Inverse Fast Fourier Transform from Fourier space to real space

        Args:
          normalise (bool): 1/V sum e^ikx so that it is same as
                            \int d3k/(2pi)^3 f(k) e^ikx
                            sum e^ikx if normalise = False
        """
        if self.mode != 'fourier-space':
            raise RuntimeError('grid is not in Fourier space')
        
        c._grid_fft_inverse(self._grid)
        if self.shifted is not None:
            self.shifted.fft_inverse(normalise=normalise)

        if normalise:
            self[:] *= 1.0/self.boxsize**3

        return self

    def save_h5(self, filename):
        """Save grid to an HDF5 file
        
        """
        import h5py
        
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
            
            if self.mode == 'real-space':
                g['fx'] = self.shifted[:]
            elif self.mode == 'fourier-space':
                g['fk'] = self.shifted[:]

        f.close()

    def dump_binary(self, filename, *, dtype='float'):
        """
        Save grid to an binary file

        Args:
          dtype (str): 'float', 'double', or 'native'
                       'native' uses Float specified at compile time
        """

        if self.mode == 'fourier-space':
            raise TypeError('dump_binary not implemented for fourier-space')

        assert(self.mode == 'real-space')

        if dtype == 'float':
            float_size = 4
        elif dtype == 'double':
            float_size = 8
        elif dtype == 'native':
            float_size = 0
        else:
            raise ValueError('Unknown dtype to write: {}'.format(dtype))

        c._grid_write_binary_real(self._grid, filename, float_size);


    def assign_density(self, xyz, *,
                       weight=None, nbar=None,
                       mas='CIC', parallel='default'):
        """Assign density to this grid
        Args:
            xyz    (array[n, 3]): position of particles in cartisian coodinate
            weight (array[n]):
            nbar:  (array[n]):    mean density without clustering
            mas:  The mass assignment scheme 'NGP', 'CIC', or 'TSC'
        """

        if self.mode != 'real-space':
            raise TypeError('grid is not in real-space mode')
    
        if self.boxsize <= 0.0:
            raise AssertionError('grid.boxsize is not set.')

        mas = mas.upper()
        
        if parallel == 'default':
            c._mass_assignment_from_array(xyz, weight, nbar,
                                          lssps._mass_assignment_scheme[mas],
                                          self._grid)
        elif parallel == 'serial' or parallel == 'atomic':
            c._mass_assignment_variations(parallel,
                                          xyz, weight, nbar,
                                          lssps._mass_assignment_scheme[mas],
                                          self._grid)
        else:
            raise TypeError('Unknown parallel = %s in Grid.assign_mass()' %
                            parallel)

        if self.shifted is not None:
            self.shifted.assign_density(xyz=xyz, weight=weight,
                                        nbar=nbar, mas=mas, parallel=parallel)


    def correct_mas(self):
        if self.mode != 'fourier-space':
            raise RuntimeError('Grid must be in fourier-space to correct_mas()')
            
        if not self.mas_corrected:
            c._mass_assignment_correct_mas(self._grid)
        else:
            warnings.warn('MAS already corrected. Neglect correct_mas()',
                          RuntimeWarning)
            return

        if (self.shifted is not None) and (self.interlaced == False):
            self.shifted.correct_mas()

        self.mas_corrected = True

    @property
    def mode(self):
        """Current mode of the grid:
        'real-space' or 'fourier-space'
        """
        return c._grid_get_mode(self._grid)

    @mode.setter
    def mode(self, m):
        if m == 'real-space':
            c._grid_set_mode(self._grid, 1)
        elif m == 'fourier-space':
            c._grid_set_mode(self._grid, 2)
        else:
            raise TypeError('grid.mode must be real-space or fourier-space')

        if self.shifted is not None:
            self.shifted.mode = m

        return m

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
        """Returns (w_sum, w2_sum, nw2_sum)"""
        ss = c._grid_get_sums(self._grid)
        return c._grid_get_sums(self._grid)

    @sums.setter
    def sums(self, values):
        c._grid_set_sums(self._grid, values[0], values[1], values[2], values[3])

    @property
    def w_sum(self):
        """ sum_i w_i"""
        ss = self.sums
        return ss[0]

    @property
    def w2_sum(self):
        """ sum_i w_i^2"""
        ss = self.sums
        return ss[1]

    @property
    def nw2_sum(self):
        """ sum_i nbar w_i^2"""
        ss = self.sums
        return ss[2]

    @property
    def n_mas(self):
        """degree of mass assignment scheme 1 for NGP, 2 for CIC, 3 for TSC"""
        return c._grid_get_nmas(self._grid)

    @n_mas.setter
    def n_mas(self, value):
        c._grid_set_nmas(self._grid, value)

    @property
    def pk_normalisation(self):
        """Normalisation factor for power spectrum"""
        return c._grid_get_pk_normalisation(self._grid)

    @pk_normalisation.setter
    def pk_normalisation(self, value):
        c._grid_set_pk_normalisation(self._grid, value)

    @property
    def shot_noise(self):
        """Value of shot noise"""
        return c._grid_get_param(self._grid, 'shot_noise')

    @shot_noise.setter
    def shot_noise(self, value):
        c._grid_set_param_double(self._grid, 'shot_noise', value)

    def shot_noise_n(self, n):
        """
        Shot noise for n>0 window multipoles
        """
        w2_sum = c._grid_get_w2(self._grid, n)
        norm = self.nw2_sum
        return w2_sum/norm
        

def empty(nc, boxsize, x0=None, offset=0.0, *, interlacing=True):
    """Return a new empty grid (uninitialised).
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

    return grid



def zeros(nc, boxsize, x0=None, offset=0.0, *, interlacing=True):
    """Return a new grid initialised with zeros.
    Args:
        nc (int): number of grids per dimension
        boxsize: length of the box on a side
        x0: corner of the box
        offset: offset of the grid points from the box corner
                give a tuple of 2 floats to set both offsets for
                the main and shifted grid
        interlacing: attach shifted grid which is shifted by half gridspacing
    """
    
    grid = empty(nc, boxsize, x0, offset, interlacing=interlacing)
    grid.clear()

    return grid

def zeros_like(grid):
    return zeros(grid.nc, grid.boxsize, grid.x0, grid.offset,
                 interlacing=(grid.shifted is not None))

def empty_like(grid):
    return empty(grid.nc, grid.boxsize, grid.x0, grid.offset,
                 interlacing=(grid.shifted is not None))

def load_h5(filename):
    """
    Args:
        filename (str): file name of an HDF file
    """
    import h5py
    
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

        grid.shifted.sums = (sums[0], sums[1], sums[2], n)
        grid.shifted.n_mas = n_mas

        a_shifted = g['fx'][:]
        c._grid_load_fx_from_array(grid.shifted._grid, a_shifted)
    
    f.close()

    return grid


def resize_fourier(grid, nc_new):
    grid_new= empty(nc_new, grid.boxsize, interlacing=False)
    
    c._grid_resize_fourier(grid._grid, grid_new._grid)

    return grid_new


def k(kind, nc, boxsize):
    grid = empty(nc, boxsize)
    grid.mode = 'fourier-space'
    
    if kind == 'mag':
        c._grid_create_kmag(grid._grid)
    else:
        if isinstance(kind, int) and 0 <= kind < 3:
            c._grid_create_k(grid._grid, kind)
        else:
            raise ValueError('kind is not 0, 1, 2, or mag')

    return grid

def mu2(nc, boxsize, *, axis=2):
    """
    Returns:
      new grid of mu2 = (k[axis]/k)^2
    """
    grid = empty(nc, boxsize, interlacing=False)

    if isinstance(axis, int) and 0 <= axis < 3:
        c._grid_set_mu2(grid._grid, axis)
    else:
        raise ValueError('axis must be an integer 0, 1, or 2: {}'.format(axis))

    return grid

def power_spectrum3D(nc, boxsize, k, P):
    """
    Create a 3D grid of P(k)

    Args:
      nc (int): number of grids per dimension
      boxsize (float): length of the grid box on a side
      k (array): array of k
      P (array): array of P(k)
    """

    if k.ndim != 1 or P.ndim != 1:
        raise TypeError('Expected a 1-dimensional array for k and P: {} {}'.format(k.shape, P.shape))

    if k.shape[0] != P.shape[0]:
        raise TypeError('Expected arrays of same size for k and P: {} {}'.format(k.shape, P.shape))
        
    
    grid = zeros(nc, boxsize, interlacing=False)
    
    c._grid_set_power3d(k, P, grid._grid)

    return grid

def dump_vector_binary(filename, grid_x, grid_y, grid_z):
    c._grid_write_vector_binary_real(filename,
                                     grid_x._grid,
                                     grid_y._grid,
                                     grid_z._grid)

