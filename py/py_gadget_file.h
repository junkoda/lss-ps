#ifndef PY_GADGET_FILE_H
#define PY_GADGET_FILE_H 1

#include "Python.h"

PyMODINIT_FUNC
py_gadget_file_module_init();

PyObject* py_gadget_file_read_np(PyObject* self, PyObject* args);
PyObject* py_gadget_file_read_header(PyObject* self, PyObject* args);

#endif
