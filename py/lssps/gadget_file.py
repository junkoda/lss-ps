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
    d['m'] = (1.0e10*m for m in mass)
    d['a'] = scale_factor
    d['redshift'] = redshift
    d['nfiles'] = nfiles
    d['boxsize'] = boxsize
    d['omega_m'] = omega_m
    d['omega_lambda'] = omega_lambda
    d['h'] = hubble_param

    return d
