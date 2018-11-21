#ifndef PY_GRID_H
#define PY_GRID_H 1

#include "Python.h"

PyMODINIT_FUNC
py_grid_module_init();

PyObject* py_grid_alloc(PyObject* self, PyObject* args);
PyObject* py_grid_fft(PyObject* self, PyObject* args);
PyObject* py_grid_fft_inverse(PyObject* self, PyObject* args);
PyObject* py_grid_nc(PyObject* self, PyObject* args);
PyObject* py_grid_get_mode(PyObject* self, PyObject* args);
PyObject* py_grid_set_mode(PyObject* self, PyObject* args);
PyObject* py_grid_get_boxsize(PyObject* self, PyObject* args);
PyObject* py_grid_set_boxsize(PyObject* self, PyObject* args);
PyObject* py_grid_get_x0(PyObject* self, PyObject* args);
PyObject* py_grid_set_x0(PyObject* self, PyObject* args);
PyObject* py_grid_get_offset(PyObject* self, PyObject* args);
PyObject* py_grid_set_offset(PyObject* self, PyObject* args);
PyObject* py_grid_get_sums(PyObject* self, PyObject* args);
PyObject* py_grid_set_sums(PyObject* self, PyObject* args);
PyObject* py_grid_get_nmas(PyObject* self, PyObject* args);
PyObject* py_grid_set_nmas(PyObject* self, PyObject* args);
PyObject* py_grid_get_pk_normalisation(PyObject* self, PyObject* args);
PyObject* py_grid_set_pk_normalisation(PyObject* self, PyObject* args);

PyObject* py_grid_fx_asarray(PyObject* self, PyObject* args);
PyObject* py_grid_fk_asarray(PyObject* self, PyObject* args);

PyObject* py_grid_load_fx_from_array(PyObject* self, PyObject* args);

PyObject* py_grid_clear(PyObject* self, PyObject* args);
PyObject* py_grid_copy(PyObject* self, PyObject* args);

PyObject* py_grid_compute_fluctuation(PyObject* self, PyObject* args);
PyObject* py_grid_compute_fluctuation_homogeneous(PyObject* self,
						  PyObject* args);

PyObject* py_grid_get_param(PyObject* self, PyObject* args);
PyObject* py_grid_set_param_double(PyObject* self, PyObject* args);

PyObject* py_grid_resize_fourier(PyObject* self, PyObject* args);

PyObject* py_grid_create_k(PyObject* self, PyObject* args);
PyObject* py_grid_create_kmag(PyObject* self, PyObject* args);
PyObject* py_grid_set_mu2(PyObject* self, PyObject* args);
PyObject* py_grid_set_power3d(PyObject* self, PyObject* args);

PyObject* py_grid_write_binary_real(PyObject* self, PyObject* args);
PyObject* py_grid_write_vector_binary_real(PyObject* self, PyObject* args);
#endif
