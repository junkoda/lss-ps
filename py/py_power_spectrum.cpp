#include "msg.h"
#include <string>
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


//static void py_power_spectrum_free(PyObject *obj);


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

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->nmodes[0]);
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

PyObject* py_power_spectrum_Pl_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_ps;
  int l;

  if(!PyArg_ParseTuple(args, "Oi", &py_ps, &l)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  const int nd=1;
  npy_intp dims[]= {ps->n};

  switch(l) {
  case 0:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p0[0]);
  case 1:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p1[0]);
  case 2:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p2[0]);
  case 3:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p3[0]);
  case 4:
    return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, &ps->p4[0]);
  }

  PyErr_SetString(PyExc_ValueError, "Unknow multipole moment");
  return NULL;
}

PyObject* py_power_spectrum_compute_plane_parallel(PyObject* self,
						   PyObject* args)
{
  // _power_spectrum_compute_multipoles(k_min, k_max, dk,
  //  grid, subtract_shotnoise, correct_mas)

  double k_min, k_max, dk;
  PyObject *py_grid;
  int subtract_shotnoise, correct_mas;
  int line_of_sight;

  if(!PyArg_ParseTuple(args, "dddOiii",
		       &k_min, &k_max, &dk,
		       &py_grid, &subtract_shotnoise, &correct_mas,
		       &line_of_sight)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  PowerSpectrum* const ps=
    multipole_compute_plane_parallel(k_min, k_max, dk,
				     grid,
				     subtract_shotnoise, correct_mas,
				     line_of_sight);
  
  return PyCapsule_New(ps, "_PowerSpectrum", py_power_spectrum_free);
}

PyObject* py_power_spectrum_compute_power_multipoles(PyObject* self,
						     PyObject* args)
{
  // _power_spectrum_compute_power_multipoles(k_min, k_max, dk,
  //  grid, subtract_shotnoise, correct_mas)

  double k_min, k_max, dk;
  PyObject *py_grid;
  int subtract_shotnoise, correct_mas;
  int line_of_sight;

  if(!PyArg_ParseTuple(args, "dddOiii",
		       &k_min, &k_max, &dk,
		       &py_grid, &subtract_shotnoise, &correct_mas,
		       &line_of_sight)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  PowerSpectrum* const ps=
    multipole_compute_power_multipoles(k_min, k_max, dk,
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

PyObject* py_power_spectrum_compute_yamamoto(PyObject* self,
					     PyObject* args)
{
  // compute Yamamoto-Scoccimarro or Yamamoto-Bianchi power spectrum
  // _power_spectrum_compute_yamamoto(k_min, k_max, dk,
  //  _grid, _grid2, _grid4, subtract_shotnoise, correct_mas)
  //
  // _grid4:
  //   Yamamoto-Scoccimarro is used if _delta4 is None
  //   Yamamoto-Bianchi otherwise
  
  double k_min, k_max, dk;
  PyObject *py_grid, *py_grid2, *py_grid4;
  int subtract_shotnoise, correct_mas;

  
  if(!PyArg_ParseTuple(args, "dddOOOii",
		       &k_min, &k_max, &dk,
		       &py_grid, &py_grid2, &py_grid4,
		       &subtract_shotnoise, &correct_mas)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Grid const * const grid2=
    (Grid const *) PyCapsule_GetPointer(py_grid2, "_Grid");
  py_assert_ptr(grid2);


  // Yamamoto-Scoccimaro is grid4 is None
  if(py_grid4 == Py_None) {
    PowerSpectrum* pss= 
      multipole_compute_yamamoto_scoccimarro(k_min, k_max, dk, 
					     grid, grid2,
					     subtract_shotnoise, correct_mas);

    return PyCapsule_New(pss, "_PowerSpectrum", py_power_spectrum_free);
  }

  // Yamamoto-Bianchi estimator
  Grid const * const grid4=
    (Grid const *) PyCapsule_GetPointer(py_grid4, "_Grid");
  py_assert_ptr(grid4);

  PowerSpectrum* psb= 
  multipole_compute_yamamoto_bianchi(k_min, k_max, dk, 
				     grid, grid2, grid4,
				     subtract_shotnoise, correct_mas);
  return PyCapsule_New(psb, "_PowerSpectrum", py_power_spectrum_free);
}

PyObject* py_power_spectrum_compute_yamamoto_odd(PyObject* self,
						 PyObject* args)
{
  // compute dipole using Yamamoto estimator
  // _power_spectrum_compute_yamamoto_1(k_min, k_max, dk,
  //  _grid, _grid2, _grid4, subtract_shotnoise, correct_mas)
  //

  PyObject *py_ps;
  double k_min, k_max, dk;
  PyObject *py_grid, *py_grid1, *py_grid3;
  int subtract_shotnoise, correct_mas;

  
  if(!PyArg_ParseTuple(args, "OdddOOOii", &py_ps,
		       &k_min, &k_max, &dk,
		       &py_grid, &py_grid1, &py_grid3, 
		       &subtract_shotnoise, &correct_mas)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Grid const * const grid1=
    (Grid const *) PyCapsule_GetPointer(py_grid1, "_Grid");
  py_assert_ptr(grid1);

  if(py_grid3 == Py_None) { // Dipole only
    PowerSpectrum* ps1=
      multipole_compute_yamamoto_odd_multipoles(k_min, k_max, dk,
						grid, grid1, 0,
						subtract_shotnoise,
						correct_mas);

    py_assert_ptr(ps->p1.size() == ps1->p1.size());
    ps->p1= ps1->p1;

    delete ps1;
  }
  else { // Dipole and tripole
    Grid const * const grid3=
      (Grid const *) PyCapsule_GetPointer(py_grid3, "_Grid");
    py_assert_ptr(grid3);

    PowerSpectrum* ps3=
      multipole_compute_yamamoto_odd_multipoles(k_min, k_max, dk,
						grid, grid1, grid3,
						subtract_shotnoise,
						correct_mas);

    py_assert_ptr(ps->p1.size() == ps3->p1.size());
    ps->p1= ps3->p1;
    ps->p3= ps3->p3;
  }

  Py_RETURN_NONE;
}

