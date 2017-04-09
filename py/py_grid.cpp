//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
#include "grid.h"
#include "error.h"
#include "py_grid.h"
#include "py_assert.h"
#include <iostream>
using namespace std;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#ifdef DOUBLEPRECISION
#define NPY_FLOAT_TYPE NPY_DOUBLE
#define NPY_COMPLEX_TYPE NPY_COMPLEX128
#else
#define NPY_FLOAT_TYPE NPY_FLOAT
#define NPY_COMPLEX_TYPE NPY_COMPLEX64
#endif

PyMODINIT_FUNC
py_grid_module_init()
{
  msg_printf(msg_debug, "Module py_grid initialised for numpy array.\n");
  import_array();

  return NULL;
}


static void py_grid_free(PyObject *obj);


PyObject* py_grid_alloc(PyObject* self, PyObject* args)
{
  // _grid_alloc(nc)
  // Create a new grid object

  int nc;
  if(!PyArg_ParseTuple(args, "i", &nc)) {
    return NULL;
  }

  //msg_printf(msg_verbose, "Allocating a grid nc = %d\n", nc);

  Grid* grid= 0;
  try {
    grid= new Grid(nc);
  }
  catch(MemoryError) {
    PyErr_SetString(PyExc_MemoryError, "Grid memory error");
    return NULL;
  }

  return PyCapsule_New(grid, "_Grid", py_grid_free);
}

void py_grid_free(PyObject *obj)
{
  // Delete the grid object; called automatically by Python
  Grid* const grid= (Grid*) PyCapsule_GetPointer(obj, "_Grid");
  py_assert_void(grid);

  delete grid;
}

PyObject* py_grid_fft(PyObject* self, PyObject* args)
{
  // Fourier transform the grid
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  try {
    grid->fft();
  }
  catch(AssertionError) {
    PyErr_SetString(PyExc_MemoryError, "Grid FFT error");
  }

  Py_RETURN_NONE;
}

PyObject* py_grid_nc(PyObject* self, PyObject* args)
{
  // _grid_nc(_cat)
  // Return the number of grid points per dimension
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("n", grid->nc);
}

PyObject* py_grid_mode(PyObject* self, PyObject* args)
{
  // _grid_mode(_cat)
  // Return the mode of the grid real-space / fourier-space
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  if(grid->mode == GridMode::real_space)
    return Py_BuildValue("s", "real-space");
  else if(grid->mode == GridMode::fourier_space)
    return Py_BuildValue("s", "fourier-space");

  Py_RETURN_NONE;
}

PyObject* py_grid_fx_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nd= 3;
  npy_intp nc= grid->nc;
  npy_intp ncz= 2*(nc/2 + 1);
  
  npy_intp dims[]= {nc, nc, nc};
  npy_intp strides[]= {(npy_intp) (sizeof(Float)*nc*ncz),
		       (npy_intp) (sizeof(Float)*ncz),
		       (npy_intp) (sizeof(Float))};
  return PyArray_New(&PyArray_Type, nd, dims, NPY_FLOAT_TYPE, strides,
  grid->fx, 0, 0, 0);
}


PyObject* py_grid_fk_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nd= 3;
  const npy_intp nc= grid->nc;
  const int nckz= nc/2 + 1;
  npy_intp dims[]= {nc, nc, nckz};

  return PyArray_SimpleNewFromData(nd, dims, NPY_COMPLEX_TYPE, grid->fk);
}

PyObject* py_grid_clear(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->clear();

  Py_RETURN_NONE;
}

PyObject* py_grid_compute_fluctuation(PyObject* self, PyObject* args)
{
  PyObject *py_grid_data, *py_grid_rand;
  if(!PyArg_ParseTuple(args, "OO", &py_grid_data, &py_grid_rand)) {
    return NULL;
  }

  Grid* const grid_data=
    (Grid*) PyCapsule_GetPointer(py_grid_data, "_Grid");
  py_assert_ptr(grid_data);

  Grid const * const grid_rand=
    (Grid const *) PyCapsule_GetPointer(py_grid_rand, "_Grid");
  py_assert_ptr(grid_rand);

  grid_compute_fluctuation(*grid_data, *grid_rand);

  Py_RETURN_NONE;
}


PyObject* py_grid_compute_fluctuation_homogeneous(PyObject* self,
						  PyObject* args)
{
  PyObject *py_grid_data, *py_grid_rand;
  if(!PyArg_ParseTuple(args, "OO", &py_grid_data, &py_grid_rand)) {
    return NULL;
  }

  Grid* const grid_data=
    (Grid*) PyCapsule_GetPointer(py_grid_data, "_Grid");
  py_assert_ptr(grid_data);


  grid_compute_fluctuation_homogeneous(*grid_data);

  Py_RETURN_NONE;
}


