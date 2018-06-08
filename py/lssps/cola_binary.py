import lssps._lssps as c
import numpy as np

def load_particles_header(filename):
    h = c._cola_binary_load_particles_header(filename)
    d = {'np': h[0],
         'boxsize': h[1],
         'm': h[2],
         'omega_m': h[3],
         'h': h[4],
         'a': h[5],
         'z': h[6]}

    return d

def load_particles(filename):
    """
    Returns:
      xv (array): [number of particles, 6]
      header (dict): header
    """
    h = load_particles_header(filename)

    a = np.empty((h['np'], 6))
    c._cola_binary_load_particles(filename, a)

    return a, h
