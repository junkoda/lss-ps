#include <iostream>
#include "msg.h"
#include "catalogue.h"
#include "grid.h"
#include "mass_assignment.h"
#include "error.h"
#include "py_mass_assignment.h"
#include "py_assert.h"
#include "py_util.h"

using namespace std;


PyObject* py_mass_assignment_from_array(PyObject* self, PyObject* args)
{
  // _mass_assignment_from_array(xyz, weight, nbar, mas, _grid)
  // Args:
  //     xyz:    The array of xyz
  //     weight: The array of weights, can be None
  //     nbar:   The array of nbar, can be None
  //     mas (int): 1 for NGP, 2 for CIC, ...
  //     _grid: a _Grid pointer

  PyObject *py_xyz, *py_weight, *py_nbar;
  PyObject *py_grid;
  int mas;
  Py_ssize_t zero[]= {(Py_ssize_t) 0};

  if(!PyArg_ParseTuple(args, "OOOiO",
		       &py_xyz, &py_weight, &py_nbar, &mas,
		       &py_grid))
    return NULL;

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  //
  // Decode array information
  //

  // xyz array
  Py_buffer xyz;
  if(PyObject_GetBuffer(py_xyz, &xyz, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;
  
  if(xyz.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array for xyz");
    return NULL;
  }

  // weight array
  Py_buffer weight;
  if(py_weight == Py_None) {
    weight.buf= nullptr;
    weight.suboffsets= nullptr;
    weight.strides= zero;
  }
  else {
    if(PyObject_GetBuffer(py_weight, &weight, PyBUF_FORMAT | PyBUF_FULL_RO)
       == -1) return NULL;

    if(weight.ndim != 1) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected a 1-dimensional array for weight");
      return NULL;
    }
  }

  // nbar
  Py_buffer nbar;
  if(py_nbar == Py_None) {
    nbar.buf= nullptr;
    nbar.suboffsets= nullptr;
    nbar.strides= zero;
  }
  else {
    if(PyObject_GetBuffer(py_nbar, &nbar, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
      return NULL;

    if(nbar.ndim != 1) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected a 1-dimensional array for nbar");
      return NULL;
    }
  }


  if(xyz.suboffsets || weight.suboffsets || nbar.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_mass_assignment_from_array cannot handle array with suboffsets");
    return NULL;
  }

  if(weight.buf && (weight.shape[0] != xyz.shape[0])) {
    PyErr_SetString(PyExc_TypeError, "Length of xyz and weight arrays differ");
    return NULL;
  }

  if(nbar.buf && (nbar.shape[0] != xyz.shape[0])) {
    PyErr_SetString(PyExc_TypeError, "Length of xyz and nbar arrays differ");
    return NULL;
  }

  if(strcmp(xyz.format, "d") == 0) {
    if(xyz.strides[1] != sizeof(double)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    mass_assignment_from_array<double>((double*) xyz.buf,
				       xyz.strides[0],
				       (double*) weight.buf,
				       weight.strides[0],
				       (double*) nbar.buf,
				       nbar.strides[0],
				       xyz.shape[0],
				       mas, true, grid);
  }
  else if(strcmp(xyz.format, "f") == 0) {
    if(xyz.strides[1] != sizeof(float)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    mass_assignment_from_array<float>((float*) xyz.buf,
				       xyz.strides[0],
				       (float*) weight.buf,
				       weight.strides[0],
				       (float*) nbar.buf,
				       nbar.strides[0],
				       xyz.shape[0],
				       mas, true, grid);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "Expected an array of floats or doubles");
    return NULL;
  }


  PyBuffer_Release(&xyz);
  if(weight.buf)
    PyBuffer_Release(&weight);
  if(nbar.buf)
    PyBuffer_Release(&nbar);


  Py_RETURN_NONE;
}

