//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
#include <iostream>
#include "power_spectrum.h"
#include "py_assert.h"

using namespace std;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#ifdef DOUBLEPRECISION
#define NPY_FLOAT_TYPE NPY_DOUBLE
#else
#define NPY_FLOAT_TYPE NPY_FLOAT
#endif

PyMODINIT_FUNC
py_power_spectrum_module_init()
{
  import_array();

  return NULL;
}


static void py_power_spectrum_free(PyObject *obj);

PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args)
{
  // _power_spectrum_alloc(kmin, kmax, dk)
  // Create a new power spectrum object
  double kmin, kmax, dk;
  if(!PyArg_ParseTuple(args, "ddd", &kmin, &kmax, &dk)) {
    return NULL;
  }

  PowerSpectrum* const ps= new PowerSpectrum(kmin, kmax, dk);

  return PyCapsule_New(ps, "_PowerSpectrum", py_power_spectrum_free);
}

void py_power_spectrum_free(PyObject *obj)
{
  // Delete the power spectrum object
  // Called automatically by Python
  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(obj, "_PowerSpectrum");
  py_assert_void(ps);

  delete ps;
}


PyObject* py_power_spectrum_len(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  cerr << "_len\n";

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  return Py_BuildValue("i", ps->n);
}


PyObject* py_power_spectrum_k_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  const int nd=1;
  npy_intp dims[]= {ps->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, ps->k_hist);
}


PyObject* py_power_spectrum_P0_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  cerr << "_P0\n";

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  const int nd=1;
  npy_intp dims[]= {ps->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, ps->P0_hist);
}

PyObject* py_power_spectrum_P2_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  const int nd=1;
  npy_intp dims[]= {ps->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, ps->P2_hist);
}

PyObject* py_power_spectrum_P4_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  const int nd=1;
  npy_intp dims[]= {ps->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, ps->P4_hist);
}

PyObject* py_power_spectrum_compute_multipoles(PyObject* self, PyObject* args)
{
  // _power_spectrum_compute_multipoles(grid, subtract_shotnoise, neff, _ps)
  PyObject *py_grid, *py_ps;
  int subtract_shotnoise;
  double neff;

  if(!PyArg_ParseTuple(args, "OidO",
		       &py_grid, &subtract_shotnoise, &neff, &py_ps)) {
    return NULL;
  }

  cerr << "ps\n";
  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  cerr << "grid\n";
  Grid* const grid= (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  power_spectrum_compute_multipoles(grid, ps, subtract_shotnoise, neff);

  Py_RETURN_NONE;  
}
