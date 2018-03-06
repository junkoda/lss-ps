#include <iostream>
#include "msg.h"
#include "catalogue.h"
#include "grid.h"
#include "mean_density.h"
#include "error.h"
#include "py_mass_assignment.h"
#include "py_assert.h"
#include "py_util.h"

using namespace std;


PyObject* py_mean_density_from_grid(PyObject* self, PyObject* args)
{
  // _mass_assignment_from_array(_grid, xyz, nbar)
  // Args:
  //     _grid: a _Grid pointer
  //     xyz:    The array of xyz
  //     nbar:   The array of nbar, can be None

  PyObject *py_xyz, *py_nbar;
  PyObject *py_grid;
  double fac;

  if(!PyArg_ParseTuple(args, "OdOO",
		       &py_grid, &fac,
		       &py_xyz, &py_nbar))
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

  // nbar
  Py_buffer nbar;
  if(PyObject_GetBuffer(py_nbar, &nbar, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;

  if(nbar.ndim != 1) {
    PyErr_SetString(PyExc_TypeError,
		    "Expected a 1-dimensional array for nbar");
    return NULL;
  }


  if(xyz.suboffsets || nbar.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_mass_assignment_from_array cannot handle array with suboffsets");
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

    mean_density_from_grid<double>(grid,
				   fac,
				   xyz.shape[0],
				   (double*) xyz.buf,
				   xyz.strides[0],
				   (double*) nbar.buf,
				   nbar.strides[0]);
  }
  else if(strcmp(xyz.format, "f") == 0) {
    if(xyz.strides[1] != sizeof(float)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    mean_density_from_grid<float>(grid,
				  fac,
				  xyz.shape[0],
				  (float*) xyz.buf,
				  xyz.strides[0],
				  (float*) nbar.buf,
				  nbar.strides[0]);
  }

  else {
    PyErr_SetString(PyExc_TypeError, "Expected an array of doubles");
    return NULL;
  }


  PyBuffer_Release(&xyz);
  PyBuffer_Release(&nbar);


  Py_RETURN_NONE;
}

