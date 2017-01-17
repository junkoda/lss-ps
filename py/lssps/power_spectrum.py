import lssps._lssps as c

class PowerSpectrum:
    def __init__(self, kmin, kmax, dk):
        self._ps = c._power_spectrum_alloc(kmin, kmax, dk)
        

    def __len__(self):
        return c._power_spectrum_len(self._ps)

    def update(self):
        self.k = c._power_spectrum_k_asarray(self._ps)
        self.P0 = c._power_spectrum_P0_asarray(self._ps)

def compute_multipoles(grid, *, kmin=0.0, kmax=1.0, dk=0.01,
                       subtract_shotnoise = True, neff= -1.6):
    ps = PowerSpectrum(kmin, kmax, dk)

    c._power_spectrum_compute_multipoles(grid._grid,
                                         int(subtract_shotnoise), neff, ps._ps)
    ps.update()
    
    return ps
