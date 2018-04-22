import lssps._lssps as c
from lssps.power_spectrum import PowerSpectrum
import lssps

class Model:
    def __init__(self, _model):
        """
        Args:
          _model: _Model pointer to C++ Model object
        """
        self._model = _model

    def __call__(self, k, mu):
        """
        Returns
          Model power spectrum P(k, mu)
        """
        return c._model_call(self._model, k, mu)

    def compute_multipoles_discrete(self, nc, boxsize, *,
                                    k_min=0.0, k_max=1.0, dk=0.01,
                                    mas_window_nc=0, mas_window_deg=0):
        """
        Compute PowerSpectrum on descrete k grid
        """
        _ps = c._model_compute_discrete_multipoles(self._model, nc, boxsize,
                                                   k_min, k_max, dk,
        		                        mas_window_nc, mas_window_deg)

        return PowerSpectrum(_ps)
    
    def apply_window_3d(self, grid_window, k):
        """
        Apply 3D window function grid to the model

        Args:
          k (float): 

        Returns:
          convoluted power spectrum
          P~(k) = \int d3k'/(2pi)^3 |W(k - k')|^2 P(k')
          for k=(k, 0, 0)
        """
        return c._model_apply_window_3d(self._model, grid_window._grid, k)

    def grid(self, nc, boxsize):
        g = lssps.grid.empty(nc, boxsize, interlacing=False)
        c._model_create_grid(self._model, g._grid)

        return g


def linear(filename, omega_m, z=0.0, *, b=1.0,
           redshift_space=True, sigma_v=0.0):

    _model= c._model_linear_alloc(filename, omega_m, z, b,
                                  redshift_space, sigma_v)

    return Model(_model)


