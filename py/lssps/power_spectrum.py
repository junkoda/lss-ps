import lssps._lssps as c
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


