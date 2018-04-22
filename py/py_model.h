#ifndef PY_MODEL_H
#define PY_MODEL_H 1

#include "Python.h"

PyObject* py_model_linear_alloc(PyObject* self, PyObject* args);
PyObject* py_model_call(PyObject* self, PyObject* args);
PyObject* py_model_compute_discrete_multipoles(PyObject* self, PyObject* args);
PyObject* py_model_apply_window_3d(PyObject* self, PyObject* args);
#endif
