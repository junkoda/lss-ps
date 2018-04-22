#include <iostream>
#include "cosmology.h"
#include "py_cosmology.h"
#include "py_assert.h"

PyObject* py_cosmology_init(PyObject* self, PyObject* args)
{
  // py_cosmology_init(omega_m, z_max, n)
  double omega_m, z_max;
  int n;
  
  if(!PyArg_ParseTuple(args, "ddi", &omega_m, &z_max, &n))
    return NULL;

  cosmology_init(omega_m, z_max, n);

  Py_RETURN_NONE;
}

PyObject* py_cosmology_distance_redshift(PyObject* self, PyObject* args)
{
  // Compute redshift from distance
  // _cosmolgy_distance_redshift(d, z)
  // Args:
  //    d: input distance
  //    z: output redshift
  PyObject *py_d, *py_z;
  if(!PyArg_ParseTuple(args, "OO", &py_d, &py_z))
    return NULL;

  Py_buffer a_d, a_z;
  if(PyObject_GetBuffer(py_d, &a_d, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    return NULL;
  if(PyObject_GetBuffer(py_z, &a_z, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    return NULL;
  
  if(a_d.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array for distance");
    return NULL;
  }

  if(a_z.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array for redshift");
    return NULL;
  }

  if(strcmp(a_d.format, "d") != 0) {
    PyErr_SetString(PyExc_TypeError,
		    "distance is expected to be an array of double");
    return NULL;
  }
  if(strcmp(a_z.format, "d") != 0) {
    PyErr_SetString(PyExc_TypeError,
		    "redshift is expected to be an array of double");
    return NULL;
  }

  py_assert_ptr(a_z.shape[0] == a_d.shape[0]);
  size_t n = a_z.shape[0];
  
  double const * const d= (double*) a_d.buf;
  double * const z= (double*) a_z.buf;

  for(size_t i=0; i<n; ++i) {
    z[i]= cosmology_distance_redshift(d[i]);
  }

  PyBuffer_Release(&a_d);
  PyBuffer_Release(&a_z);

  Py_RETURN_NONE;
}


PyObject* py_cosmology_compute_comoving_distance(PyObject* self, PyObject* args)
{
  double z;
  if(!PyArg_ParseTuple(args, "d", &z))
    return NULL;

  double a= 1.0/(1.0 + z);
  double d= cosmology_compute_comoving_distance(a);

  return Py_BuildValue("d", d);
}

//
// growth factor
//
PyObject* py_cosmology_growth_D(PyObject* self, PyObject* args)
{
  // _growth_D(a, omega_m)
  double a, omega_m;
  if(!PyArg_ParseTuple(args, "dd", &a, &omega_m)) {
    return NULL;
  }

  return Py_BuildValue("d", cosmology_growth_D(a, omega_m));
}

PyObject* py_cosmology_growth_rate_f(PyObject* self, PyObject* args)
{
  // _growth_f(a, omega_m)
  double a, omega_m;
  if(!PyArg_ParseTuple(args, "dd", &a, &omega_m)) {
    return NULL;
  }

  return Py_BuildValue("d", cosmology_growth_f(a, omega_m));
}
