#ifndef PY_CATALOGUE_H
#define PY_CATALOGUE_H 1

#include "Python.h"

PyMODINIT_FUNC
py_catalogue_module_init();

PyObject* py_catalogue_alloc(PyObject* self, PyObject* args);
PyObject* py_catalogue_len(PyObject* self, PyObject* args);
PyObject* py_catalogue_read_text(PyObject* self, PyObject* args);

#endif
