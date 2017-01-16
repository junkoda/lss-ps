//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
#include "catalogue.h"
#include "grid.h"
#include "mass_assignment.h"
#include "error.h"
#include "py_mass_assignment.h"
#include "py_assert.h"

PyObject* py_mass_assignment_cic(PyObject* self, PyObject* args)
{
  // _mass_assignment_cic(_cat, x0, boxsize, _grid)
  PyObject *py_cat, *py_grid, *py_x0;
  double boxsize;

  if(!PyArg_ParseTuple(args, "OOdO", &py_cat, &py_x0, &boxsize, &py_grid)) {
    return NULL;
  }

  Catalogue const * const cat=
    (Catalogue const *) PyCapsule_GetPointer(py_cat, "_Catalogue");
  py_assert_ptr(cat);

  Grid* const grid= (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Float x0[]= {0, 0, 0};

  mass_assignment_cic(cat, x0, boxsize, grid);

  Py_RETURN_NONE;
}

