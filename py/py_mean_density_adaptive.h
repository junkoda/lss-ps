#ifndef PY_MEAN_DENSITY_ADAPTIVE_H
#define PY_MEAN_DENSITY_ADAPTIVE_H 1

#include "Python.h"

//PyMODINIT_FUNC
//py_mean_density_adaptive_module_init();


PyObject* py_kdpoints_alloc(PyObject* self, PyObject* args);
PyObject* py_kdpoints_len(PyObject* self, PyObject* args);
PyObject* py_kdpoints_from_array(PyObject* self, PyObject* args);
PyObject* py_kdpoints_density_as_array(PyObject* self, PyObject* args);

PyObject* py_kdtree_alloc(PyObject* self, PyObject* args);

PyObject* py_kdtree_compute_rk(PyObject* self, PyObject* args);
PyObject* py_kdtree_update_node_statistics(PyObject* self, PyObject* args);

PyObject* py_mean_density_adaptive_estimate(PyObject* self, PyObject* args);
PyObject* py_mean_density_average_estimate(PyObject* self, PyObject* args);


#endif
