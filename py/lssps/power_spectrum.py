import lssps._lssps as c

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
        self._ps = _ps
        self.k = c._power_spectrum_k_asarray(self._ps)
        self.P0 = c._power_spectrum_P0_asarray(self._ps)
        self.P2 = c._power_spectrum_P2_asarray(self._ps)
        self.P4 = c._power_spectrum_P4_asarray(self._ps)
        
        self.n = c._power_spectrum_len(self._ps)

    def __len__(self):
        """Number of bins"""
        return self.n




