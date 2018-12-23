"""
Module for computing power spectrum.

PowerSpectrum: class for the estimated power spectrum

compute_plane_parallel(grid_delta, *,
                       k_min=0.0, k_max=1.0, dk=0.01,
                       subtract_shotnoise=True,
                       correct_mas= True, line_of_sight=2)

def compute_yamamoto(grid_delta, kind, *,
                     k_min=0.0, k_max=1.0, dk=0.01,
                     subtract_shotnoise=True,
                     correct_mas=True,
                     grid_moment=None, grid2=None, grid4=None)
kind = 'scoccimarro' or 'bianchi'
"""

import lssps._lssps as c
import lssps.yamamoto
import numbers
import numpy as np

class PowerSpectrum:
    """
    Power spectrum computed in bins
    Attributes:
        k (numpy.array): mean wave number in bins [h/Mpc]
        P0 (numpy.array): P0 monopole [(1/h Mpc)^3]
        P2 (numpy.array): P2 quadrupole
        P4 (numpy.array): P4 hexadecapole
        nmodes (numpy.array): number of independent k modes in bins
    """
    def __init__(self, _ps):
        """
        Arg:
            _ps: _PowerSpectrum pointer
        """
        
        self._ps = _ps
        self.k = c._power_spectrum_k_asarray(self._ps)
        self.nmodes = c._power_spectrum_nmodes_asarray(self._ps).astype(float)
        
        self.P0 = c._power_spectrum_P0_asarray(self._ps)
        self.P2 = c._power_spectrum_P2_asarray(self._ps)
        self.P4 = c._power_spectrum_P4_asarray(self._ps)

        self.P1 = c._power_spectrum_Pl_asarray(self._ps, 1)
        self.P3 = c._power_spectrum_Pl_asarray(self._ps, 3)

        self.n = c._power_spectrum_len(self._ps)
        self.shot_noise = c._power_spectrum_shotnoise(self._ps)

    def __len__(self):
        """Number of bins"""
        return self.n

    def __getitem__(self, index):
        if isinstance(index, slice):
            k = self.k[index]
            n = len(k)
            
            a = np.zeros(n, 5)
            a[:, 0] = k
            a[:, 1] = self.P0[index]
            a[:, 2] = self.P2[index]
            a[:, 3] = self.P4[index]
            a[:, 4] = self.nmodes[index]

            return a

        elif isinstance(index, numbers.Integral):
            i=index
            return (self.k[i], self.nmodes[i],
                    self.P0[i], self.P2[i], self.P4[i])

        else:
            raise TypeError('Index must be an integer')


def compute_plane_parallel(grid_delta, *,
                           k_min=0.0, k_max=1.0, dk=0.01,
                           subtract_shotnoise=True,
                           correct_mas= True, line_of_sight=2):
    """
    Args:
        delta_grid (Grid): Grid object for delta(k)
        k_min (float): lower bound of k binning [h/Mpc]
        k_max (float): upper bound of k binning [h/Mpc]
        dk (float):    bin width [h/Mpc]
        subtract_shotnoise (bool)
        correct_mas (bool)
        line_of_sight (int): line of sight direction for multipoles
                             0,1,2 for x,y,z, respectively
    Returns:
        PowerSpectrum object
    """
    
    if grid_delta.mode != 'fourier-space':
        grid_delta.fft()

    if grid_delta.shifted is not None:
        if grid_delta.interlaced == False:
            grid_delta.interlace()

    correct_mas = correct_mas and (not grid_delta.mas_corrected)

    _ps = c._power_spectrum_compute_plane_parallel(k_min, k_max, dk,
                              grid_delta._grid,
                              subtract_shotnoise, correct_mas,
                              line_of_sight)

    return PowerSpectrum(_ps)

    
compute_yamamoto = lssps.yamamoto.compute_yamamoto

def compute_power_multipoles(grid, *,
                             k_min=0.0, k_max=1.0, dk=0.01,
                             subtract_shotnoise=False,
                             correct_mas= False, line_of_sight=2):
    """
    Compute multipoles of 3D power spectrum P(kvec)
    Args:
        grid (Grid):   Grid object for P(k)
        k_min (float): lower bound of k binning [h/Mpc]
        k_max (float): upper bound of k binning [h/Mpc]
        dk (float):    bin width [h/Mpc]
        subtract_shotnoise (bool)
        correct_mas (bool)
        line_of_sight (int): line of sight direction for multipoles
                             0,1,2 for x,y,z, respectively
    Returns:
        PowerSpectrum object
    """
    
    if grid.mode != 'fourier-space':
        raise ValueError('grid.mode must be fourier-space: %s' % grid.mode)

    correct_mas = correct_mas and (not grid.mas_corrected)

    _ps = c._power_spectrum_compute_power_multipoles(k_min, k_max, dk,
                              grid._grid,
                              subtract_shotnoise, correct_mas,
                              line_of_sight)

    return PowerSpectrum(_ps)


def compute_discrete_multipoles(grid, *,
                                k_min=0.0, k_max=1.0, dk=0.01,
                                subtract_shotnoise=True,
                                correct_mas= True, line_of_sight=2):
    """
    Args:
        grid (Grid): Grid object for delta(k)
        k_min (float): lower bound of k binning [h/Mpc]
        k_max (float): upper bound of k binning [h/Mpc]
        dk (float):    bin width [h/Mpc]
        subtract_shotnoise (bool)
        correct_mas (bool)
        line_of_sight (int): line of sight direction for multipoles
                             0,1,2 for x,y,z, respectively
    Returns:
        PowerSpectrum object
    """

    correct_mas = correct_mas and (not grid.mas_corrected)

    if grid.mode != 'fourier-space':
        grid.fft()

    if grid.shifted is not None:
        if grid.interlaced == False:
            grid.interlace()

    _ps = c._power_spectrum_compute_discrete_multipoles(1,
                              k_min, k_max, dk,
                              grid._grid,
                              subtract_shotnoise, correct_mas,
                              line_of_sight)
        
    return PowerSpectrum(_ps)


def compute_discrete_power_multipoles(grid, *,
                                k_min=0.0, k_max=1.0, dk=0.01,
                                subtract_shotnoise=False,
                                correct_mas=False, line_of_sight=2):
    """
    Args:
        grid (Grid): Grid object of delta(k) or P(k)
        kind (str): The content of grid 'delta' or 'power' 
        k_min (float): lower bound of k binning [h/Mpc]
        k_max (float): upper bound of k binning [h/Mpc]
        dk (float):    bin width [h/Mpc]
        subtract_shotnoise (bool)
        correct_mas (bool)
        line_of_sight (int): line of sight direction for multipoles
                             0,1,2 for x,y,z, respectively
    Returns:
        PowerSpectrum object
    """

    correct_mas = correct_mas and (not grid.mas_corrected)

    if grid.mode != 'fourier-space':
        raise ValueError('grid.mode must be fourier-space: %s' % grid.mode)

    _ps = c._power_spectrum_compute_discrete_multipoles(0,
                              k_min, k_max, dk,
                              grid._grid,
                              subtract_shotnoise, correct_mas,
                              line_of_sight)
        
    return PowerSpectrum(_ps)


def compute_2d_average(grid, *, k_min=0.0, k_max=1.0, dk=0.01, nmu=10):
    """
    compute 2D average (k, mu) of the 3D grid
    """
    nk= round((k_max - k_min)/dk)

    if nk <= 0:
        raise ValueError('No k bin in k range %e %e with dk=%e' %
                         (k_min, k_max, dk))

    if nmu <= 0:
        raise ValueError('nmu must be positive: %d' % nmu)

    nbin= nk*nmu;
    k = np.empty(nbin)
    mu = np.empty_like(k)
    P = np.empty_like(k)
    nmodes = np.empty_like(k)
    
    c._power_spectrum_compute_2d_average(k_min, dk, nk, nmu, grid._grid,
                                         nmodes, k, mu, P)

    k = k.reshape((nk, nmu))
    mu = mu.reshape((nk, nmu))
    P = P.reshape((nk, nmu))
    nmodes = nmodes.reshape((nk, nmu))

    d = {}
    d['k'] = k
    d['mu'] = mu
    d['P'] = P
    d['nmodes'] = nmodes

    return d
    
def compute_power_2d(grid, *, k_min=0.0, k_max=1.0, dk=0.01, nmu=10,
                     line_of_sight=2,
                     subtract_shotnoise=True, correct_mas=True):
    """
    Compute 2D (auto) power spectrum from the grid

    Args:
      grid (Grid): Grid objects in Fourier space (TODO)
      k_min (float): minimum edge of k binning
      k_max (float): number of bin will be round((k_max - k_min)/dk)
      dk (float): k bin width
      nmu (int): number of mu bins
      line_of_sight (int): i = 0,1, or, 2 for mu = k[i]/|k|
      shot_noise (bool or float): subtract the given shot_noise if it is float
                                  subtract grid.shot_noise if True
      correct_mas (bool): correct for the mass assignment window function

    Returns:
      d (dict):
        d['nmodes']: nmodes[ik, imu] number of independent modes
        d['k']:      mean k in the 2D bin
        d['mu']:     mean mu in the 2D  bin
        d['P']:      power spectrum in the 2D bin
    """
    
    if isinstance(subtract_shotnoise, bool):
        if subtract_shotnoise == True:
            shot_noise = grid.shot_noise
        else:
            shot_noise = 0.0
    elif isinstance(subtract_shotnoise, float):
        shot_noise = subtract_shotnoise
    else:
        raise TypeError('Unrecognised type for subtract_shotnoise; '
                        'it must be a boolean or the shot noise value in float')

    return compute_cross_power_2d(grid, grid, 're',
                                  k_min=k_min, k_max=k_max, dk=dk, nmu=nmu,
                                  shot_noise=shot_noise,
                                  correct_mas=correct_mas)


def compute_cross_power_2d(grid1, grid2, kind='re', *,
                           k_min=0.0, k_max=1.0, dk=0.01, nmu=10,
                           line_of_sight=2,
                           shot_noise=0.0, correct_mas=True):
    """
    Compute 2D cross-power spectrum from grid1 and grid2
    Re(delta1 delta2^*) or Im(delta1 delta2^*)

    Args:
      grid1, grid2 (Grid): Grid objects in Fourier space (TODO)
      kind (str): 're' or 'im', Real or Imaginary part of grid1 grid2^*
      k_min (float): minimum edge of k binning
      k_max (float): number of bin will be round((k_max - k_min)/dk)
      dk (float): k bin width
      nmu (int): number of mu bins
      line_of_sight (int): i = 0,1, or, 2 for mu = k[i]/|k|
      shot_noise (float): subtract this as a shot noise
      correct_mas (bool): correct for the mass assignment window function

    Returns:
      d (dict):
        d['nmodes']: nmodes[ik, imu] number of independent modes
        d['k']:      mean k in the 2D bin
        d['mu']:     mean mu in the 2D  bin
        d['P']:      power spectrum in the 2D bin
    
      k, mu, P are 0 if no modes are in the bin
    """

    assert(nmu > 0)
    
    if kind == 're':
        real_imag = 0
    elif kind == 'im':
        real_imag = 1
    else:
        raise ValueError("kind must be 're' or 'im': %s" % kind)
    
    nk = round((k_max - k_min)/dk)
    if nk <= 0:
        raise ValueError('No k bins for the given k_min, k_max, dk')

    # FFT if necessary
    if grid1.mode != 'fourier-space':
        grid1.fft()
    if grid2.mode != 'fourier-space':
        grid2.fft()

    
    
    n = nk*nmu
    
    nmodes = np.zeros(n)
    k = np.zeros(n)
    mu = np.zeros(n)
    ps2d = np.zeros(n)

    c._power2d_compute(grid1._grid, grid2._grid, real_imag,
                       k_min, dk, nk, nmu,
		       shot_noise, correct_mas, line_of_sight,
		       nmodes, k, mu, ps2d)
    
    d = {}
    d['nmodes'] = nmodes.reshape((nk, nmu))
    d['k'] = k.reshape((nk, nmu))
    d['mu'] = mu.reshape((nk, nmu))
    d['P'] = ps2d.reshape((nk, nmu))

    return d
