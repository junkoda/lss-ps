#include "grid.h"
#include "py_assert.h"
#include "lognormal.h"
#include "py_lognormal.h"

PyObject* py_lognormal_convert_P_to_delta_k(PyObject* self, PyObject* args)
{
  // _lognormal_convert_P_to_delta_k(_grid, seed, fix_amplitude)
  PyObject *py_grid;
  unsigned long seed;
  int fix_amplitude;
  
  if(!PyArg_ParseTuple(args, "Oki", &py_grid, &seed, &fix_amplitude)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  py_assert_ptr(grid->mode == GridMode::fourier_space);

  convert_P_to_delta_k(grid, seed, fix_amplitude);

  Py_RETURN_NONE;
}
