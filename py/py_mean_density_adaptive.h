#ifndef PY_MEAN_DENSITY_ADAPTIVE_H
#define PY_MEAN_DENSITY_ADAPTIVE_H 1

#include "Python.h"

PyObject* py_kdpoints_alloc(PyObject* self, PyObject* args);
PyObject* py_kdpoints_from_array(PyObject* self, PyObject* args);
PyObject* py_mean_density_adaptive_estimate(PyObject* self, PyObject* args);
#endif
