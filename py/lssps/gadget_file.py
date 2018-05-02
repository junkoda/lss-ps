import math
import numpy as np
import lssps._lssps as c

def read_np(filename):
    return c._gadget_file_read_np(filename)

def read_header(filename):
    """
    Read header in Gadget binary file
    
    Args:
      filename (str): file name
    
    Returns:
      dictionary of parameters
      np: number of particles in this file
      np_total: number of total particles in the simulation
      nfiles: number of files that one snapshot is split into
      m: mass of a particle [1/h Solar mass]
      a: scale factor of the snapshot
      redshift: redshift of the snapshot
      omega_m: Omega matter at present
      omega_lambda: Omega Lambda at present
      h: Hubble constant H0 in units of 100 km/s/Mpc
    """
    
    np, np_total, mass, scale_factor, redshift, nfiles, boxsize, omega_m, omega_lambda, hubble_param = c._gadget_file_read_header(filename)
    d = {}

    d['np'] = np
    d['np_total'] = np_total
    d['m'] = mass
    d['a'] = scale_factor
    d['redshift'] = redshift
    d['nfiles'] = nfiles
    d['boxsize'] = boxsize
    d['omega_m'] = omega_m
    d['omega_lambda'] = omega_lambda
    d['h'] = hubble_param

    return d

class GadgetFile:
    def __init__(self, filename):
        self._gf = c._gadget_file_alloc(filename)
        self.header = read_header(filename)
        self.rescale_velocity = True

    def __enter__(self):
        c._gadget_file_open(self._gf)
        return self


    def __exit__(self, exc_type, exc_value, traceback):
        c._gadget_file_close(self._gf)

    def __getitem__(self, components):
        ncol_x = 0
        ncol_v = 0

        for comp in components:
            if comp == 'x':
                ncol_x = 3
            elif comp == 'v':
                ncol_v = 3
            else:
                raise ValueError('Unknown component for GadgetFile; must be x or v: %s' % comp)

        ncol = ncol_x + ncol_v
        
        assert(self.header['np'][0] == 0)
        nrow = self.header['np'][1]

        a = np.empty((nrow, ncol), dtype='f')
        n = self.header['np'][1]

        if ncol_x:
            c._gadget_file_read(self._gf, 'x', 0, n, a[:, 0:3])
            icol=3
        if ncol_v:
            c._gadget_file_read(self._gf, 'v', 0, n, a[:, ncol_x:(ncol_x + 3)])
            if self.rescale_velocity:
                vfac = math.sqrt(self.header['a'])
                a[:, ncol_x:(ncol_x + 3)] *= vfac

        return a
        
    
def open(filename, *, rescale_velocity=True):
    gf = GadgetFile(filename)
    gf.rescale_velocity = rescale_velocity
    return gf
