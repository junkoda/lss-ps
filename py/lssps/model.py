import lssps._lssps as c
from lssps.power_spectrum import PowerSpectrum

class Model:
    def __init__(self, _model):
        self._model = _model

    def __call__(self, k, mu):
        return c._model_call(self._model, k, mu)

    def compute_multipoles_discrete(self, nc, boxsize, *,
                                    k_min=0.0, k_max=1.0, dk=0.01,
                                    mas_window_nc=0, mas_window_deg=0):
        """
        Compute PowerSpectrum on descrete k grid
        """
        print('DEBUG mas_window_deg', mas_window_deg)
        _ps = c._model_compute_discrete_multipoles(self._model, nc, boxsize,
                                                   k_min, k_max, dk,
        		                           mas_window_nc, mas_window_deg)

        return PowerSpectrum(_ps)
    
    def apply_mas():
        pass
    
    def apply_window3d(grid_window):
        pass


def linear(filename, omega_m, z=0.0, *, b=1.0,
           redshift_space=True, sigma_v=0.0):

    _model= c._model_linear_alloc(filename, omega_m, z, b,
                                  redshift_space, sigma_v)

    return Model(_model)


