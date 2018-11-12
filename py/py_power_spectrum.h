#ifndef PY_POWER_SPECTRUM_H
#define PY_POWER_SPECTRUM_H 1

#include "Python.h"

PyMODINIT_FUNC
py_power_spectrum_module_init();

void py_power_spectrum_free(PyObject *obj);

//PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_len(PyObject* self, PyObject* args);

PyObject* py_power_spectrum_k_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_nmodes_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P0_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P2_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P4_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_Pl_asarray(PyObject* self, PyObject* args);

//PyObject* py_power_spectrum_compute_multipoles(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_compute_plane_parallel(PyObject* self,
						   PyObject* args);
PyObject* py_power_spectrum_compute_yamamoto(PyObject* self, PyObject* args);

PyObject* py_power_spectrum_compute_yamamoto_odd(PyObject* self,
						 PyObject* args);

PyObject* py_power_spectrum_shotnoise(PyObject* self, PyObject* args);

PyObject* py_power_spectrum_compute_power_multipoles(PyObject* self,
						     PyObject* args);



#endif
