#ifndef PY_GADGET_FILE_H
#define PY_GADGET_FILE_H 1

#include "Python.h"

PyObject* py_gadget_file_read_np(PyObject* self, PyObject* args);
PyObject* py_gadget_file_read_header(PyObject* self, PyObject* args);

PyObject* py_gadget_file_alloc(PyObject* self, PyObject* arg);
PyObject* py_gadget_file_open(PyObject* self, PyObject* args);
PyObject* py_gadget_file_close(PyObject* self, PyObject* args);
PyObject* py_gadget_file_read(PyObject* self, PyObject* args);

#endif
