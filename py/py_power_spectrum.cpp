#include "msg.h"
#include <iostream>
#include "power_spectrum.h"
#include "grid.h"
#include "multipole.h"
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


/*
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
*/

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

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->k[0]);
}

PyObject* py_power_spectrum_nmodes_asarray(PyObject* self, PyObject* args)
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_INT, &ps->nmodes[0]);
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

  //for(int i=0; i<ps->n; ++i) {
    //fprintf(stderr, "%e %e 0\n", ps->k[i], (*ps->p0[i]);
    //}

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p0[0]);
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p2[0]);
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

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p4[0]);
}

PyObject* py_power_spectrum_compute_plane_parallel(PyObject* self,
						   PyObject* args)
{
  // _power_spectrum_compute_multipoles(k_min, k_max, dk, nmu,
  //  grid, subtract_shotnoise, correct_mas)
  
  double k_min, k_max, dk;
  int nmu;
  PyObject *py_grid;
  int subtract_shotnoise, correct_mas;
  int line_of_sight;
  
  if(!PyArg_ParseTuple(args, "dddiOiii",
		       &k_min, &k_max, &dk, &nmu,
		       &py_grid, &subtract_shotnoise, &correct_mas,
		       &line_of_sight)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  PowerSpectrum* const ps=
    multipole_compute_plane_parallel(k_min, k_max, dk, nmu,
				     grid,
				     subtract_shotnoise, correct_mas,
				     line_of_sight);
  
  return PyCapsule_New(ps, "_PowerSpectrum", py_power_spectrum_free);
}

PyObject* py_power_spectrum_shotnoise(PyObject* self, PyObject* args)
{
  PyObject *py_ps;

  if(!PyArg_ParseTuple(args, "O", &py_ps)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  return Py_BuildValue("d", ps->shot_noise);
}
