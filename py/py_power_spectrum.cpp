//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
//#include "catalogue.h"
//#include "error.h"
//#include "py_catalogue.h"
#include "power_spectrum.h"
#include "py_assert.h"

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

PyObject* py_double_array1(double* const p, const int n)
{
  const int nd=1;
  npy_intp dims[]= {n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE, p);
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE, ps->k_hist);
}

PyObject* py_power_spectrum_P0_asarray(PyObject* self, PyObject* args)
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE, ps->P0_hist);
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE, ps->P2_hist);
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE, ps->P4_hist);
}
