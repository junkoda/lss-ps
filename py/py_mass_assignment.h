#ifndef PY_MASS_ASSIGNMENT
#define PY_MASS_ASSIGNMENT 1

#include "Python.h"

PyObject* py_mass_assignment(PyObject* self, PyObject* args);
PyObject* py_mass_assignment_from_array(PyObject* self, PyObject* args);

#endif
