import lssps
import lssps._lssps as c
from numbers import Number
import warnings
import h5py
impo

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
        grid.interlacing: 'done' if interlace() has been called. 
                           Reset to None by clear()
    """

    def __init__(self, nc, boxsize, x0=None, offset=None):
        # _grid is the pointer to a C++ Grid object
        self._grid = c._grid_alloc(nc)
        self.shifted = None
        self.boxsize = boxsize
        self.interlaced = False

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
        self.interlaced = False
        self.mas_corrected = False
        
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
        self.interlacing = True

        return self

    def fft(self):
        """Fast Fourier Transform from real space to Fourier space"""
        if self.mode != 'real-space':
            raise RuntimeError('grid is not in real space')

        c._grid_fft(self._grid)
        if self.shifted is not None:
            self.shifted.fft()

        return self

    def fft_inverse(self):
        """Inverse Fast Fourier Transform from Fourier space to real space"""
        if self.mode != 'fourier-space':
            raise RuntimeError('grid is not in Fourier space')
        
        c._grid_fft_inverse(self._grid)
        if self.shifted is not None:
            self.shifted.fft_inverse()

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
            
            if self.mode == 'real-space':
                g['fx'] = self.shifted[:]
            elif self.mode == 'fourier-space':
                g['fk'] = self.shifted[:]

        f.close()

    def assign_density(self, *,
                       cat=None,
                       xyz=None, weight=None, nbar=None,
                       mas='CIC', parallel='default'):
        """Assign density to this grid
        Args:
            mas:  The mass assignment scheme 'NGP', 'CIC', or 'TSC'
            xyz    (array[n, 3]): position of particles in cartisian coodinate
            weight (array[n]):
            nbar:  (array[n]):    mean density without clustering
        """

        if self.mode != 'real-space':
            raise TypeError('grid is not in real-space mode')
    
        if self.boxsize <= 0.0:
            raise AssertionError('grid.boxsize is not set.')

        mas = mas.upper()
        
        if xyz is not None:
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

                                              
        else:
            RuntimeError('xyz not provided')
            
        if self.shifted is not None:
            self.shifted.assign_density(cat=cat, xyz=xyz, weight=weight,
                                        nbar=nbar, mas=mas, parallel=parallel)


    def correct_mas(self):
        if self.mode != 'fourier-space':
            raise RuntimeError('Grid must be in fourier-space to correct_mas()')
            
        if not self.mas_corrected:
            c._mass_assignment_correct_mas(self._grid)
        else:
            warnings.warn('MAS already corrected. Neglect correct_mas()',
                          RuntimeWarning)

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
        #print('setting pk_nromalisataion', value)
        c._grid_set_pk_normalisation(self._grid, value)

    @property
    def shot_noise(self):
        """Value of shot noise"""
        return c._grid_get_param(self._grid, 'shot_noise')

    @shot_noise.setter
    def shot_noise(self, value):
        c._grid_set_param_double(self._grid, 'shot_noise', value)




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

        grid.shifted.sums = (sums[0], sums[1], sums[2], n)
        grid.shifted.n_mas = n_mas

        a_shifted = g['fx'][:]
        c._grid_load_fx_from_array(grid.shifted._grid, a_shifted)
    
    f.close()

    return grid


