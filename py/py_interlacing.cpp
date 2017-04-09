#include "py_interlacing.h"


PyObject* py_interlacing(PyObject* self, PyObject* args)
{
  // _interlacing(_grid, _grid_shifted)
  PyObject *py_grid1, *py_grid2;

  if(!PyArg_ParseTuple(args, "OO", &py_grid1, &py_grid2)) {
    return NULL;
  }

  Grid* const grid1= (Grid*) PyCapsule_GetPointer(py_grid1, "_Grid");
  py_assert_ptr(grid1);

  Grid* const grid2= (Grid*) PyCapsule_GetPointer(py_grid2, "_Grid");
  py_assert_ptr(grid2);


  interlacing(grid1, grid2);

  Py_RETURN_NONE;  
}
