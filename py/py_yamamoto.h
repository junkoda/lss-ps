#ifndef PY_YAMAMOTO_H
#define PY_YAMAMOTO_H

#include "Python.h"

PyObject* py_yamamoto_compute_moment_x(PyObject* self, PyObject* args);
PyObject* py_yamamoto_compute_moment_k(PyObject* self, PyObject* args);

#endif
