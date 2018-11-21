#include "py_lognormal.h"

PyObject* py_lognormal_convert_P_to_delta_k(PyObject* self, PyObject* args)
{
  // _lognormal_convert_P_to_delta_k(_grid)
  PyObject *py_grid;
  unsigned long seed;
  int fix_ampletude;
  
  if(!PyArg_ParseTuple(args, "Oki", &py_grid, &seed, &fix_amplitude)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  py_assert_ptr(grid->mode == fourier_space);

  convert_P_to_delta_k(grid, seed, fix_amplitude);

  Py_RETURN_NONE;
}
