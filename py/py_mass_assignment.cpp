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

using namespace std;
/*
PyObject* py_mass_assignment(PyObject* self, PyObject* args)
{
  // _mass_assignment_cic(_cat, x0, boxsize, mas, _grid)
  PyObject *py_cat, *py_grid, *py_x0;
  int mas;
  double boxsize;

  if(!PyArg_ParseTuple(args, "OOdiO",
		       &py_cat, &py_x0, &boxsize, &mas, &py_grid))
    return NULL;

  Catalogue const * const cat=
    (Catalogue const *) PyCapsule_GetPointer(py_cat, "_Catalogue");
  py_assert_ptr(cat);

  Grid* const grid= (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Float x0[]= {0, 0, 0};

  mass_assignment(cat, x0, boxsize, mas, grid);

  Py_RETURN_NONE;
}
*/

PyObject* py_mass_assignment(PyObject* self, PyObject* args)
{
  PyObject *py_f, *py_grid, *py_x0;
  int mas;
  double boxsize;

  if(!PyArg_ParseTuple(args, "OOdiO",
		       &py_f, &py_x0, &boxsize, &mas, &py_grid))
    return NULL;

  CatalogueFile* const f=
    (CatalogueFile*) PyCapsule_GetPointer(py_f, "_CatalogueFile");
  py_assert_ptr(f);

  // TODO; with interlacing, py_grid will be a pair of grids
  Grid* const grid= (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  Grid* const grid_shifted= 0;
  py_assert_ptr(grid);

  // TODO set x0
  Float x0[]= {0, 0, 0};
  Float x0_shifted[]= {0, 0, 0};

  static size_t buffer_size = 100000;
  
  f->open();

  #pragma omp parallel
  {
    vector<Particle> cat;
    cat.reserve(buffer_size);

    while(true) {
      f->read(buffer_size, cat);
      if(cat.empty())
	break;
      
      mass_assignment(cat, x0, boxsize, mas, false, grid);

      if(grid_shifted)
	mass_assignment(cat, x0_shifted, boxsize, mas, false, grid_shifted);
    }
  }

  f->close();

  Py_RETURN_NONE;
}


