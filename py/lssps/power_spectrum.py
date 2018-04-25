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

        P2D (numpy.array): 2D power spectrum P(k, mu)
        nomodes2D (numpy.array): number of k modes in 2D bins
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
            a[:, 1] = self.nmodes[index]
            a[:, 1] = self.P0[index]
            a[:, 2] = self.P2[index]
            a[:, 3] = self.P4[index]

            return a

        elif isinstance(index, numbers.Integral):
            i=index
            return (self.k[i], self.nmodes[i],
                    self.P0[i], self.P2[i], self.P4[i])

        else:
            raise TypeError('Index must be an integer')


def compute_plane_parallel(grid_delta, *,
                           k_min=0.0, k_max=1.0, dk=0.01, nmu=0,
                           subtract_shotnoise=True,
                           correct_mas= True, line_of_sight=2):
    """
    Args:
        delta_grid (Grid): Grid object for delta(k)
        k_min (float): lower bound of k binning [h/Mpc]
        k_max (float): upper bound of k binning [h/Mpc]
        dk (float):    bin width [h/Mpc]
        nmu (int):   number of mu = kz/k bins for 2D power spectrum
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
