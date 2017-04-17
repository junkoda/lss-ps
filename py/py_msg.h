#ifndef MY_MSG_H
#define MY_MSG_H 1

#include "Python.h"

PyObject* py_msg_set_loglevel(PyObject* self, PyObject* args);
PyObject* py_msg_get_loglevel(PyObject* self, PyObject* args);

#endif
