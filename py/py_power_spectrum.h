#ifndef PY_POWER_SPECTRUM_H
#define PY_POWER_SPECTRUM_H

#include "Python.h"

PyMODINIT_FUNC
py_power_spectrum_module_init();


PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_len(PyObject* self, PyObject* args);

PyObject* py_power_spectrum_k_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P0_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P2_asarray(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_P4_asarray(PyObject* self, PyObject* args);

PyObject* py_power_spectrum_compute_multipoles(PyObject* self, PyObject* args);
#endif
