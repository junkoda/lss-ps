#ifndef PY_RR_H
#define PY_RR_H 1

#include "Python.h"

PyMODINIT_FUNC
py_rr_module_init();

PyObject* py_rr_compute_multipoles(PyObject* self, PyObject* args);
PyObject* py_rr_asarray(PyObject* self, PyObject* args);


#endif
