#include <iostream> //debug

#include "rr.h"
#include "py_rr.h"
#include "py_assert.h"


#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;

PyMODINIT_FUNC
py_rr_module_init()
{
  import_array();

  return NULL;
}

static void py_rr_free(PyObject *obj);



void py_rr_free(PyObject *obj)
{
  // Delete the RRMultipoles object
  // Called automatically by Python
  RRMultipoles* const rr=
    (RRMultipoles*) PyCapsule_GetPointer(obj, "_RRMultipoles");
  py_assert_void(rr);

  delete rr;
}

PyObject* py_rr_compute_multipoles(PyObject* self, PyObject* args)
{
  PyObject *py_xyz, *py_weight, *py_nbar;
  double r_max, dr;

  if(!PyArg_ParseTuple(args, "OOOdd",
		       &py_xyz, &py_weight, &py_nbar, &r_max, &dr))
    return NULL;

  RRMultipoles* const rr= new RRMultipoles(r_max, dr);

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

  Py_ssize_t zero[]= {(Py_ssize_t) 0};

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

    if(strcmp(weight.format, xyz.format) != 0) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected same type for xyz and weight");
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

    if(strcmp(nbar.format, xyz.format) != 0) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected same type for xyz and nbar");
    }
  }


  if(xyz.suboffsets || weight.suboffsets || nbar.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_rr_compute_multipoles cannot handle array with suboffsets");
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

    rr_multipoles<double>((double*) xyz.buf,
			  xyz.strides[0],
			  (double*) weight.buf,
			  weight.strides[0],
			  (double*) nbar.buf,
			  nbar.strides[0],
			  xyz.shape[0],
			  rr);
  }
  else if(strcmp(xyz.format, "f") == 0) {
    if(xyz.strides[1] != sizeof(float)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    rr_multipoles<float>((float*) xyz.buf,
			  xyz.strides[0],
			  (float*) weight.buf,
			  weight.strides[0],
			  (float*) nbar.buf,
			  nbar.strides[0],
			  xyz.shape[0],
			  rr);
  }

  return PyCapsule_New(rr, "_RRMultipoles", py_rr_free);
}

PyObject* py_rr_asarray(PyObject* self, PyObject* args)
{
  //
  // _rr_asarray(rr, l)
  // rr: _RRMultipoles pointer
  // l (int): multipole index 0 -- 4
  //
  PyObject *py_rr;
  int l;

  if(!PyArg_ParseTuple(args, "Oi", &py_rr, &l)) {
    return NULL;
  }

  RRMultipoles* const rr=
    (RRMultipoles*) PyCapsule_GetPointer(py_rr, "_RRMultipoles");
  py_assert_ptr(rr);

  const int nd=1;
  npy_intp dims[]= {rr->n};

  switch(l) {
  case 0:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, rr->rr0.data());
  case 1:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, rr->rr1.data());
  case 2:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, rr->rr2.data());
  case 3:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, rr->rr3.data());
  case 4:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, rr->rr4.data());

  }
  
  PyErr_SetString(PyExc_ValueError, "l is not 0,1,2,3,4");
  return NULL;
}
